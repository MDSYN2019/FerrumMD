use nalgebra::Vector3;
use rand::Rng;
use rand_distr::{Distribution, Normal};
use sang_md::cell_subdivision::SimulationBox;
use sang_md::lennard_jones_simulations::{
    apply_thermostat_berendsen_particles, compute_forces_particles, compute_temperature_particles,
    pbc_update, LJParameters, Particle,
};
use sang_md::molecule::io::{write_gro, write_xtc};
use sang_md::molecule::martini;
use sang_md::molecule::molecule::System;

fn create_tip3p_water_box(
    n_side: usize,
    box_length: f64,
    temperature: f64,
    mass: f64,
    n_particles: usize,
) -> Result<Vec<System>, String> {
    let spacing = box_length / n_side as f64; // spacing between molecules
    let sigma_v = (temperature / mass).sqrt();
    let normal = Normal::new(0.0, sigma_v)
        .map_err(|e| format!("failed to build normal distribution: {e}"))?;
    let mut rng = rand::rng();
    // create storage vector to store the
    let mut systems = Vec::with_capacity(n_particles);

    // Start modifying the molecule coordinates to be able to fit in the simulation box

    let mut tip3p_forcefield_itp_instance =
        martini::ItpForceField::read_itp("ff/Charmm27.ff/charmm27.ff/tip3p.itp")?;
    let shared_atom_types = martini::ItpForceField::read_atomtypes_from_itp(
        "ff/Charmm27.ff/charmm27.ff/ffnonbonded.itp",
    )?;
    tip3p_forcefield_itp_instance
        .atom_types
        .extend(shared_atom_types);

    // set up the initial oxygen and the h1/h2 positions
    let mut o = Vector3::new(0.0, 0.0, 0.0); // oxygen position
    let mut h1 = Vector3::new(0.9572, 0.0, 0.0); // h1 position
    let mut h2 = Vector3::new(-0.23999, 0.92663, 0.0); // h2 position
    let mut o_h1_h2_vec: [Vector3<f64>; 3] = [o, h1, h2];
    // initiate the coordinates
    //let mut tip3p_instance_initiated = tip3p_forcefield_itp_instance.to_system(o_h1_h2_vec);

    for ix in 0..n_side {
        for iy in 0..n_side {
            for iz in 0..n_side {
                let jitter = 0.05 * spacing;
                // this will be the center cartesian point to seed the water molecule
                let origin = Vector3::new(
                    (ix as f64 + 0.5) * spacing + rng.random_range(-jitter..jitter),
                    (iy as f64 + 0.5) * spacing + rng.random_range(-jitter..jitter),
                    (iz as f64 + 0.5) * spacing + rng.random_range(-jitter..jitter),
                );

                // copy an instance of the coordinate initiated tip3p molecule, modifying the position with the position coordinate above
                let mut o_h1_h2_vec: [Vector3<f64>; 3] = [o + origin, h1 + origin, h2 + origin];

                // initiate the coordinates
                let mut tip3p_instance_initiated = tip3p_forcefield_itp_instance
                    .to_system(&o_h1_h2_vec)
                    .expect("the tip3p build should succeed");

                // push the instance
                systems.push(tip3p_instance_initiated);
                // generate the water molecule
            }
        }
    }

    Ok(systems)
}

fn snapshot(particles: &[Particle]) -> Vec<Particle> {
    particles.to_vec()
}

fn main() {
    // Martini-style CG water-bead NVT setup (single-bead solvent model).
    let n_side = 6;
    let target_temperature = 300.0;
    let mass = 72.0;
    let sigma = 0.47;
    let epsilon = 0.2;
    let box_length = 8.0;
    let dt = 0.002;
    let nsteps = 200;
    let thermostat_tau = 0.05;

    //let mut particles =
    //    create_martini_water_box(n_side, box_length, target_temperature, mass, sigma, epsilon)?;
    //
    //let simulation_box = SimulationBox {
    //    x_dimension: box_length,
    //    y_dimension: box_length,
    //    z_dimension: box_length,
    //};
    //
    //let mut subcells = simulation_box.create_subcells(10);
    //simulation_box.store_atoms_in_cells_particles(&mut particles, &mut subcells, 10);
    //compute_forces_particles(&mut particles, box_length, &mut subcells);
    //
    //let mut frames = Vec::with_capacity((nsteps / 20) as usize + 1);
    //frames.push(snapshot(&particles));
    //
    //for step in 0..nsteps {
    //    let a_old: Vec<_> = particles.iter().map(|p| p.force / p.mass).collect();
    //
    //    for (p, a) in particles.iter_mut().zip(a_old.iter()) {
    //        p.velocity += 0.5 * a * dt;
    //        p.position += p.velocity * dt;
    //    }
    //
    //    pbc_update(&mut particles, box_length);
    //
    //    let mut subcells = simulation_box.create_subcells(10);
    //    simulation_box.store_atoms_in_cells_particles(&mut particles, &mut subcells, 10);
    //    compute_forces_particles(&mut particles, box_length, &mut subcells);
    //
    //    for p in &mut particles {
    //        let a_new = p.force / p.mass;
    //        p.velocity += 0.5 * a_new * dt;
    //    }
    //
    //    apply_thermostat_berendsen_particles(
    //        &mut particles,
    //        target_temperature,
    //        thermostat_tau,
    //        dt,
    //    );
    //
    //    if step % 100 == 0 {
    //        let temp =
    //            compute_temperature_particles(&particles, 3 * particles.len().saturating_sub(1));
    //        log::info!("step={step:4} T={temp:.2}");
    //    }
    //
    //    if step % 20 == 0 {
    //        frames.push(snapshot(&particles));
    //    }
    //}
    //
    //write_gro(
    //    "martini_water_box.gro",
    //    &particles,
    //    Vector3::new(box_length, box_length, box_length),
    //    "Martini CG water box (NVT)",
    //)?;
    //
    //write_xtc(
    //    "martini_water_box.xtc",
    //    &frames,
    //    Vector3::new(box_length, box_length, box_length),
    //    dt as f32,
    //)?;
    //
    //println!(
    //    "Wrote martini_water_box.gro and martini_water_box.xtc for {} particles and {} frames",
    //    particles.len(),
    //    frames.len()
    //);
    //
    //Ok(())
}
