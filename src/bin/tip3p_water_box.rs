use nalgebra::Vector3;
use rand::Rng;
use sang_md::lennard_jones_simulations::{self, InitOutput};
use sang_md::molecule::io::write_gro_systems;
use sang_md::molecule::martini;
use sang_md::molecule::molecule::System;

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
) {
    for step in 0..max_steps {
        for sys in systems.iter_mut() {
            lennard_jones_simulations::compute_all_bonded_forces_system(sys, box_length);
        }
        lennard_jones_simulations::compute_intermolecular_forces_systems(systems, box_length);
        let _ =
            lennard_jones_simulations::compute_electrostatic_forces_systems(systems, box_length);

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

        if max_force < force_tolerance {
            println!(
                "Minimization converged in {} steps (max |F| = {:.6})",
                step + 1,
                max_force
            );
            return;
        }
    }

    println!("Minimization reached max steps without full convergence (max_steps={max_steps})");
}

fn main() -> Result<(), String> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let n_side = 6;
    let box_length = 18.0;
    let dt = 0.001;
    let nsteps = 20;
    let minimization_steps = 200;
    let minimization_step_size = 0.0005;
    let minimization_force_tolerance = 1e-3;

    let mut systems = create_tip3p_water_box(n_side, box_length)?;
    minimize_systems(
        &mut systems,
        box_length,
        minimization_steps,
        minimization_step_size,
        minimization_force_tolerance,
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

    lennard_jones_simulations::run_md_nve_systems(
        &mut systems,
        nsteps,
        dt,
        box_length,
        "berendsen",
    );

    write_gro_systems(
        "tip3p_water_box.gro",
        &systems,
        Vector3::new(box_length, box_length, box_length),
        "TIP3P water box",
    )?;

    //write_xtc(
    //    "martini_water_box.xtc",
    //    &frames,
    //    Vector3::new(box_length, box_length, box_length),
    //    dt as f32,
    //)?;

    let atom_count: usize = systems.iter().map(|system| system.atoms.len()).sum();
    println!(
        "Wrote tip3p_water_box.gro for {} molecules ({} atoms)",
        systems.len(),
        atom_count
    );

    Ok(())
}
