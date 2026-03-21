use nalgebra::Vector3;
use rand::Rng;
use sang_md::lennard_jones_simulations::{self, InitOutput, Particle};
use sang_md::molecule::io::write_gro;
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
    let h1 = Vector3::new(0.9572, 0.0, 0.0);
    let h2 = Vector3::new(-0.23999, 0.92663, 0.0);

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

fn flatten_systems_to_particles(systems: &[System]) -> Vec<Particle> {
    systems
        .iter()
        .flat_map(|system| system.atoms.iter().cloned())
        .collect()
}

fn main() -> Result<(), String> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let n_side = 6;
    let box_length = 18.0;
    let dt = 0.001;
    let nsteps = 20;

    let mut systems = create_tip3p_water_box(n_side, box_length)?;

    let mut init_state = InitOutput::Systems(systems.clone());
    lennard_jones_simulations::set_molecular_positions_and_velocities(&mut init_state, 300.0);
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

    let particles = flatten_systems_to_particles(&systems);
    write_gro(
        "tip3p_water_box.gro",
        &particles,
        Vector3::new(box_length, box_length, box_length),
        "TIP3P water box",
    )?;

    println!(
        "Wrote tip3p_water_box.gro for {} molecules ({} atoms)",
        systems.len(),
        particles.len()
    );

    Ok(())
}
