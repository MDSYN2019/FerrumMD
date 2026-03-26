use nalgebra::Vector3;
use rand::Rng;
use sang_md::lennard_jones_simulations::{
    self, ConstraintMode, ConstraintOptions, InitOutput, SystemSimulationConfig,
};
use sang_md::molecule::io::{systems_to_particles_frame, write_gro_systems, write_xtc};
use sang_md::molecule::martini;
use sang_md::molecule::molecule::System;
use sang_md::molecule::shake_rattle::shake_rattle;
use sang_md::PmeConfig;

fn create_tip3p_water_box(n_side: usize, box_length: f64) -> Result<Vec<System>, String> {
    let spacing = box_length / n_side as f64;
    let mut rng = rand::rng();
    let mut systems = Vec::with_capacity(n_side * n_side * n_side);

    let mut tip3p_forcefield_itp_instance =
        martini::ItpForceField::read_itp("ff/Charmm27.ff/charmm27.ff/tip3p.itp")?;
    let shared_atom_types = martini::ItpForceField::read_atomtypes_from_itp(
        "ff/Charmm27.ff/charmm27.ff/ffnonbonded.itp",
    )?;
    tip3p_forcefield_itp_instance
        .atom_types
        .extend(shared_atom_types);

    let o = Vector3::new(0.0, 0.0, 0.0);
    let h1 = Vector3::new(0.09572, 0.0, 0.0);
    let h2 = Vector3::new(-0.023999, 0.092663, 0.0);

    for ix in 0..n_side {
        for iy in 0..n_side {
            for iz in 0..n_side {
                let jitter = 0.05 * spacing;
                let origin = Vector3::new(
                    (ix as f64 + 0.5) * spacing + rng.random_range(-jitter..jitter),
                    (iy as f64 + 0.5) * spacing + rng.random_range(-jitter..jitter),
                    (iz as f64 + 0.5) * spacing + rng.random_range(-jitter..jitter),
                );

                let coordinates = [o + origin, h1 + origin, h2 + origin];
                let tip3p_instance_initiated =
                    tip3p_forcefield_itp_instance.to_system(&coordinates)?;
                systems.push(tip3p_instance_initiated);
            }
        }
    }

    Ok(systems)
}

fn minimize_systems(
    systems: &mut [System],
    box_length: f64,
    max_steps: usize,
    step_size: f64,
    force_tolerance: f64,
    lj_cutoff: f64,
    pme: PmeConfig,
) {
    log::info!(
        "Starting minimization: max_steps={}, step_size={}, force_tolerance={}, lj_cutoff={}, pme(alpha={}, real_cutoff={}, kmax={})",
        max_steps,
        step_size,
        force_tolerance,
        lj_cutoff,
        pme.alpha,
        pme.real_cutoff,
        pme.kmax
    );

    for step in 0..max_steps {
        for sys in systems.iter_mut() {
            lennard_jones_simulations::compute_all_bonded_forces_system(sys, box_length);
        }
        lennard_jones_simulations::compute_intermolecular_forces_systems_cutoff(
            systems, box_length, lj_cutoff,
        );
        let _ = lennard_jones_simulations::compute_electrostatic_forces_systems_with_config(
            systems, box_length, &pme,
        );

        let mut max_force = 0.0;
        for sys in systems.iter_mut() {
            for atom in sys.atoms.iter_mut() {
                let force_norm = atom.force.norm();
                if force_norm > max_force {
                    max_force = force_norm;
                }
                atom.position += (step_size / atom.mass) * atom.force;
            }
            lennard_jones_simulations::pbc_update(&mut sys.atoms, box_length);
        }

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

    let n_side = 6;
    let target_number_density = 33.3679; // molecules / nm^3 (~1 g/cm^3 water)
    let n_molecules = (n_side * n_side * n_side) as f64;
    let box_length = (n_molecules / target_number_density).cbrt();
    let dt = 0.001;
    let nsteps = 1000;
    let trajectory_stride = 50;
    let minimization_steps = 50;
    let minimization_step_size = 0.0005;
    let minimization_force_tolerance = 5e-3;
    let minimization_lj_cutoff = (0.5 * box_length).min(1.2);
    let minimization_pme = PmeConfig {
        alpha: 3.0,
        real_cutoff: minimization_lj_cutoff,
        kmax: 4,
    };

    let mut systems = create_tip3p_water_box(n_side, box_length)?;
    log::info!(
        "TIP3P short run config: minimization_steps={}, production_steps={}, dt={}",
        minimization_steps,
        nsteps,
        dt
    );
    minimize_systems(
        &mut systems,
        box_length,
        minimization_steps,
        minimization_step_size,
        minimization_force_tolerance,
        minimization_lj_cutoff,
        minimization_pme,
    );

    let mut init_state = InitOutput::Systems(systems.clone());
    lennard_jones_simulations::set_molecular_positions_and_velocities(
        &mut init_state,
        300.0,
        box_length,
    );
    if let InitOutput::Systems(randomized) = init_state {
        systems = randomized;
    }

    let constraints_by_system: Result<Vec<_>, _> = systems
        .iter()
        .map(shake_rattle::tip3p_constraints_from_system)
        .collect();
    let constraint_options = ConstraintOptions {
        mode: ConstraintMode::SettlePreferred,
        constraints_by_system: constraints_by_system?,
        tolerance: 1e-10,
        max_iter: 100,
    };

    let cutoff = (0.5 * box_length).min(1.2);
    let run_config = SystemSimulationConfig {
        cutoff,
        neighbor_skin: 0.2,
        neighbor_rebuild_interval: 10,
        pme: PmeConfig {
            alpha: 3.5,
            real_cutoff: cutoff,
            kmax: 6,
        },
    };

    let mut frames = Vec::with_capacity((nsteps as usize / trajectory_stride) + 2);
    frames.push(systems_to_particles_frame(&systems));

    for step in 0..nsteps {
        lennard_jones_simulations::run_md_nve_systems_with_constraints_and_config(
            &mut systems,
            1,
            dt,
            box_length,
            "none",
            Some(&constraint_options),
            run_config,
        );

        if (step + 1) % trajectory_stride == 0 || step + 1 == nsteps {
            frames.push(systems_to_particles_frame(&systems));
        }
    }

    write_gro_systems(
        "tip3p_water_box.gro",
        &systems,
        Vector3::new(box_length, box_length, box_length),
        "TIP3P water box",
    )?;

    write_xtc(
        "tip3p_water_box.xtc",
        &frames,
        Vector3::new(box_length, box_length, box_length),
        (dt * trajectory_stride as f64) as f32,
    )?;

    let atom_count: usize = systems.iter().map(|system| system.atoms.len()).sum();
    println!(
        "Wrote tip3p_water_box.gro and tip3p_water_box.xtc for {} molecules ({} atoms), {} frames",
        systems.len(),
        atom_count,
        frames.len()
    );

    Ok(())
}
