use nalgebra::Vector3;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rand_distr::{Distribution, Normal};
use sang_md::cell_subdivision::SimulationBox;
use sang_md::lennard_jones_simulations::{
    apply_thermostat_berendsen_particles, compute_forces_particles, compute_temperature_particles,
    pbc_update, LJParameters, Particle,
};
use sang_md::molecule::io::{write_gro, write_xtc};
use std::env;

fn create_martini_water_box(
    n_side: usize,
    box_length: f64,
    temperature: f64,
    mass: f64,
    sigma: f64,
    epsilon: f64,
    rng: &mut StdRng,
) -> Result<Vec<Particle>, String> {
    let n_particles = n_side * n_side * n_side;
    let spacing = box_length / n_side as f64;
    let sigma_v = (temperature / mass).sqrt();
    let normal = Normal::new(0.0, sigma_v)
        .map_err(|e| format!("failed to build normal distribution: {e}"))?;

    let mut particles = Vec::with_capacity(n_particles);

    for ix in 0..n_side {
        for iy in 0..n_side {
            for iz in 0..n_side {
                let jitter = 0.05 * spacing;
                let position = Vector3::new(
                    (ix as f64 + 0.5) * spacing + rng.random_range(-jitter..jitter),
                    (iy as f64 + 0.5) * spacing + rng.random_range(-jitter..jitter),
                    (iz as f64 + 0.5) * spacing + rng.random_range(-jitter..jitter),
                );

                let velocity =
                    Vector3::new(normal.sample(rng), normal.sample(rng), normal.sample(rng));

                particles.push(Particle {
                    id: particles.len(),
                    position,
                    velocity,
                    force: Vector3::zeros(),
                    lj_parameters: LJParameters {
                        epsilon,
                        sigma,
                        number_of_atoms: 1,
                    },
                    mass,
                    energy: 0.0,
                    atom_type: 0.0,
                    charge: 0.0,
                });
            }
        }
    }

    Ok(particles)
}

fn snapshot(particles: &[Particle]) -> Vec<Particle> {
    particles.to_vec()
}

fn remove_center_of_mass_drift(particles: &mut [Particle]) {
    let total_mass: f64 = particles.iter().map(|p| p.mass).sum();
    if total_mass <= 0.0 {
        return;
    }

    let com_velocity = particles
        .iter()
        .fold(Vector3::zeros(), |acc, p| acc + p.velocity * p.mass)
        / total_mass;

    for particle in particles.iter_mut() {
        particle.velocity -= com_velocity;
    }
}

fn rescale_temperature(particles: &mut [Particle], target_temperature: f64) {
    let dof = 3 * particles.len().saturating_sub(1);
    if dof == 0 {
        return;
    }

    let current_temperature = compute_temperature_particles(particles, dof);
    if current_temperature <= f64::EPSILON {
        return;
    }

    let scale = (target_temperature / current_temperature).sqrt();
    for particle in particles.iter_mut() {
        particle.velocity *= scale;
    }
}

fn minimize_particles(
    particles: &mut Vec<Particle>,
    simulation_box: &SimulationBox,
    box_length: f64,
    cell_subdivisions: usize,
    max_steps: usize,
    step_size: f64,
    force_tolerance: f64,
) {
    let max_displacement = 0.002;
    let max_force_for_update = 5.0e3;
    log::info!(
        "Starting minimization: max_steps={}, step_size={}, force_tolerance={}",
        max_steps,
        step_size,
        force_tolerance
    );

    for step in 0..max_steps {
        let mut subcells = simulation_box.create_subcells(cell_subdivisions);
        simulation_box.store_atoms_in_cells_particles(particles, &mut subcells, cell_subdivisions);
        compute_forces_particles(particles, box_length, &mut subcells);

        let mut max_force = 0.0;
        for p in particles.iter_mut() {
            let force_norm = p.force.norm();
            if force_norm > max_force {
                max_force = force_norm;
            }

            let force_scale = if force_norm > max_force_for_update {
                max_force_for_update / force_norm
            } else {
                1.0
            };
            let mut displacement = step_size * p.force * force_scale;
            let displacement_norm = displacement.norm();
            if displacement_norm > max_displacement {
                displacement *= max_displacement / displacement_norm;
            }
            p.position += displacement;
        }

        pbc_update(particles, box_length);

        if step == 0 || (step + 1) % 25 == 0 || step + 1 == max_steps {
            log::info!(
                "Minimization step {:>4} | max |F| = {:.6}",
                step + 1,
                max_force
            );
        }
        if max_force < force_tolerance {
            log::info!(
                "Minimization converged in {} steps (max |F| = {:.6})",
                step + 1,
                max_force
            );
            return;
        }
    }

    log::warn!("Minimization reached max steps without full convergence (max_steps={max_steps})");
}

fn main() -> Result<(), String> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    // Martini-style CG water-bead NVT setup (single-bead solvent model).
    let n_side = 6;
    let target_temperature = 300.0;
    let mass = 72.0;
    let sigma = 0.47;
    let epsilon = 0.2;
    let box_length = 8.0;
    // Use a conservative equilibration step, then a larger production step.
    let dt_equil = 0.0002;
    let dt_prod = 0.0005;
    let nsteps_equil = 500;
    let nsteps_prod = 1500;
    let equil_thermostat_tau = 0.02;
    let prod_thermostat_tau = 0.1;
    let minimization_steps = 100;
    let minimization_step_size = 1.0e-5;
    let minimization_force_tolerance = 5e-2;
    let cell_subdivisions = 10;
    let total_steps = nsteps_equil + nsteps_prod;

    let seed = env::var("MARTINI_SEED")
        .ok()
        .and_then(|raw| raw.parse::<u64>().ok())
        .unwrap_or(20260417);
    let mut rng = StdRng::seed_from_u64(seed);
    log::info!("martini_water_box RNG seed={seed}");

    let mut particles = create_martini_water_box(
        n_side,
        box_length,
        target_temperature,
        mass,
        sigma,
        epsilon,
        &mut rng,
    )?;
    remove_center_of_mass_drift(&mut particles);
    rescale_temperature(&mut particles, target_temperature);

    let simulation_box = SimulationBox {
        x_dimension: box_length,
        y_dimension: box_length,
        z_dimension: box_length,
    };
    minimize_particles(
        &mut particles,
        &simulation_box,
        box_length,
        cell_subdivisions,
        minimization_steps,
        minimization_step_size,
        minimization_force_tolerance,
    );
    for p in particles.iter_mut() {
        p.velocity = Vector3::zeros();
    }
    rescale_temperature(&mut particles, target_temperature);

    let mut subcells = simulation_box.create_subcells(cell_subdivisions);
    simulation_box.store_atoms_in_cells_particles(&mut particles, &mut subcells, cell_subdivisions);
    compute_forces_particles(&mut particles, box_length, &mut subcells);

    let mut frames = Vec::with_capacity((total_steps / 20) as usize + 1);
    frames.push(snapshot(&particles));

    for step in 0..total_steps {
        let dt = if step < nsteps_equil {
            dt_equil
        } else {
            dt_prod
        };
        let a_old: Vec<_> = particles.iter().map(|p| p.force / p.mass).collect();

        for (p, a) in particles.iter_mut().zip(a_old.iter()) {
            p.velocity += 0.5 * a * dt;
            p.position += p.velocity * dt;
        }

        pbc_update(&mut particles, box_length);

        let mut subcells = simulation_box.create_subcells(cell_subdivisions);
        simulation_box.store_atoms_in_cells_particles(
            &mut particles,
            &mut subcells,
            cell_subdivisions,
        );
        compute_forces_particles(&mut particles, box_length, &mut subcells);

        for p in &mut particles {
            let a_new = p.force / p.mass;
            p.velocity += 0.5 * a_new * dt;
        }

        let thermostat_tau = if step < nsteps_equil {
            equil_thermostat_tau
        } else {
            prod_thermostat_tau
        };
        apply_thermostat_berendsen_particles(
            &mut particles,
            target_temperature,
            thermostat_tau,
            dt,
        );
        if step % 50 == 0 {
            remove_center_of_mass_drift(&mut particles);
        }

        if step % 100 == 0 {
            let temp =
                compute_temperature_particles(&particles, 3 * particles.len().saturating_sub(1));
            let phase = if step < nsteps_equil { "equil" } else { "prod" };
            log::info!("step={step:4} phase={phase} T={temp:.2}");
        }

        if step % 20 == 0 {
            frames.push(snapshot(&particles));
        }
    }

    write_gro(
        "martini_water_box.gro",
        &particles,
        Vector3::new(box_length, box_length, box_length),
        "Martini CG water box (NVT)",
    )?;

    write_xtc(
        "martini_water_box.xtc",
        &frames,
        Vector3::new(box_length, box_length, box_length),
        dt_equil as f32,
    )?;

    println!(
        "Wrote martini_water_box.gro and martini_water_box.xtc for {} particles and {} frames",
        particles.len(),
        frames.len()
    );

    Ok(())
}
