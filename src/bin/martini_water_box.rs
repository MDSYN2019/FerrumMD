use nalgebra::Vector3;
use rand::Rng;
use rand_distr::{Distribution, Normal};
use sang_md::cell_subdivision::SimulationBox;
use sang_md::lennard_jones_simulations::{
    apply_thermostat_berendsen_particles, compute_forces_particles, compute_temperature_particles,
    pbc_update, LJParameters, Particle,
};
use sang_md::molecule::io::{write_gro, write_xtc};

fn create_martini_water_box(
    n_side: usize,
    box_length: f64,
    temperature: f64,
    mass: f64,
    sigma: f64,
    epsilon: f64,
) -> Result<Vec<Particle>, String> {
    let n_particles = n_side * n_side * n_side;
    let spacing = box_length / n_side as f64;
    let sigma_v = (temperature / mass).sqrt();
    let normal = Normal::new(0.0, sigma_v)
        .map_err(|e| format!("failed to build normal distribution: {e}"))?;

    let mut rng = rand::rng();
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

                let velocity = Vector3::new(
                    normal.sample(&mut rng),
                    normal.sample(&mut rng),
                    normal.sample(&mut rng),
                );

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

fn main() -> Result<(), String> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    // Martini-style CG water-bead NVT setup (single-bead solvent model).
    let n_side = 6;
    let target_temperature = 300.0;
    let mass = 72.0;
    let sigma = 0.47;
    let epsilon = 0.2;
    let box_length = 8.0;
    // A conservative CG integration step helps avoid occasional LJ blow-ups
    // from rare close contacts in this simple demo setup.
    let dt = 0.005;
    let nsteps = 2000;
    let thermostat_tau = 0.02;

    let mut particles =
        create_martini_water_box(n_side, box_length, target_temperature, mass, sigma, epsilon)?;
    remove_center_of_mass_drift(&mut particles);
    rescale_temperature(&mut particles, target_temperature);

    let simulation_box = SimulationBox {
        x_dimension: box_length,
        y_dimension: box_length,
        z_dimension: box_length,
    };

    let mut subcells = simulation_box.create_subcells(10);
    simulation_box.store_atoms_in_cells_particles(&mut particles, &mut subcells, 10);
    compute_forces_particles(&mut particles, box_length, &mut subcells);

    let mut frames = Vec::with_capacity((nsteps / 20) as usize + 1);
    frames.push(snapshot(&particles));

    for step in 0..nsteps {
        let a_old: Vec<_> = particles.iter().map(|p| p.force / p.mass).collect();

        for (p, a) in particles.iter_mut().zip(a_old.iter()) {
            p.velocity += 0.5 * a * dt;
            p.position += p.velocity * dt;
        }

        pbc_update(&mut particles, box_length);

        let mut subcells = simulation_box.create_subcells(10);
        simulation_box.store_atoms_in_cells_particles(&mut particles, &mut subcells, 10);
        compute_forces_particles(&mut particles, box_length, &mut subcells);

        for p in &mut particles {
            let a_new = p.force / p.mass;
            p.velocity += 0.5 * a_new * dt;
        }

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
            log::info!("step={step:4} T={temp:.2}");
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
        dt as f32,
    )?;

    println!(
        "Wrote martini_water_box.gro and martini_water_box.xtc for {} particles and {} frames",
        particles.len(),
        frames.len()
    );

    Ok(())
}
