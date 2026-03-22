/*

=========================================================
 Molecular Dynamics Simulation Framework (Rust)
 Based on:
 "A First Encounter with the Hartree-Fock Self-Consistent Field Method"
 https://apt.scitation.org/doi/abs/10.1119/10.0002644?journalCode=ajp
=========================================================

🔧 Particle Model
-----------------
Each `Particle` struct contains:
- Position: Vector3<f64>
- Velocity: Vector3<f64>
- Force:    Vector3<f64>
- Mass:     f64
- Lennard-Jones Parameters: { sigma, epsilon }

📦 Initialization
-----------------
- Positions initialized randomly inside a cubic simulation box.
- Velocities initialized via the Maxwell-Boltzmann distribution:
    > At thermal equilibrium, particle velocities follow MB statistics.
    > This ensures that kinetic energy corresponds to the target temperature.

🌡️ Temperature Calculation
---------------------------
- Temperature derived from kinetic energy via equipartition theorem:
    T = (2/3) * (KE / N)
  Where:
    KE = total kinetic energy
    N  = number of particles

🌍 Periodic Boundary Conditions (PBC)
-------------------------------------
- Implemented as: xi ← xi mod L
- Mimics an infinite system by wrapping particles across boundaries.
- Prevents artificial wall effects and confinement artifacts.

💥 Lennard-Jones Potential
---------------------------
- Pairwise interaction calculated using:
    V(r) = 4ε [ (σ/r)^12 – (σ/r)^6 ]
- Applied Lorentz–Berthelot mixing rules for ε and σ between different species.

🧊 Velocity Rescaling Thermostat
--------------------------------
- Thermostat applied via:
    λ = sqrt(T_target / T_current)
- Rescales all velocities to control system temperature.
- Simple and stable; not strictly canonical (NVT) but effective for equilibration.

=========================================================

📌 Next Implementation Steps
----------------------------
- [ ] Force calculation from potential (∇V)
- [ ] Center-of-mass velocity removal
- [x] Minimum Image Convention (implemented)
- [x] Energy conservation check (NVE tested)
      → still to analyze for NVT/NPT cases
- [ ] Advanced thermostats:
      → Berendsen (smooth control) - Somewhat done
      → Langevin (stochastic, ensemble-correct) -
- [ ] Radial Distribution Function (RDF)

*/
extern crate assert_type_eq;

// The compiler will look for the module's code in the following places:

// src/error.rs
// src/molcule.rs
// src/parameters.rs
pub mod cell;
pub mod error;
pub mod molecule;
pub mod parameters;
#[cfg(feature = "python")]
mod python;
#[path = "quantum/quantum_chem.rs"]
pub mod quantum_chemistry;
pub mod thermostat_barostat;

use std::collections::HashSet;

// Use when importing the finished minimization modulexo
//use sang_md::lennard_jones_simulations::{self, compute_total_energy_and_print};

#[inline]
pub fn lennard_jones_force_scalar(r: f64, sigma: f64, epsilon: f64) -> f64 {
    // F(r) magnitude along r-hat; positive = repulsive
    // d/dr 4ε[(σ/r)^12 - (σ/r)^6]  =>  24ε [2(σ^12/r^13) - (σ^6/r^7)]
    if r <= 0.0 {
        return 0.0;
    }
    let sr = sigma / r;
    let sr2 = sr * sr;
    let sr6 = sr2 * sr2 * sr2;
    let sr12 = sr6 * sr6;
    24.0 * epsilon * (2.0 * sr12 - sr6) / r
}

#[inline]
fn coulomb_prefactor() -> f64 {
    // Reduced-unit Coulomb prefactor.
    // In SI units this would be 1/(4*pi*epsilon_0).
    1.0
}

#[derive(Copy, Clone, Debug)]
pub struct PmeConfig {
    pub alpha: f64,
    pub real_cutoff: f64,
    pub kmax: i32,
}

impl Default for PmeConfig {
    fn default() -> Self {
        Self {
            alpha: 0.35,
            real_cutoff: 9.0,
            kmax: 4,
        }
    }
}

#[inline]
fn safe_norm(x: f64) -> f64 {
    if x < 1e-12 {
        1e-12
    } else {
        x
    }
}

#[inline]
fn dedup_permutation(v: &mut Vec<Vec<i32>>) {
    let mut seen: HashSet<Vec<i32>> = HashSet::new();

    v.retain(|inner| {
        let mut key = inner.clone();
        key.sort_unstable(); // canonicalize: [2,1] -> [1,2]
        seen.insert(key) // true if first time
    });
}

pub mod cell_subdivision {
    use crate::cell::cell::CellList;
    use crate::lennard_jones_simulations::Particle;
    use crate::molecule::molecule::System;
    use nalgebra::Vector3;
    /*
    Create the subcells required for efficiently computing the intermolecular interactions
    within a certain radius cutoff rather than working with all the atoms in the system
     */

    fn wrap_index(i: isize, n: usize) -> usize {
        let n = n as isize;
        (((i % n) + n) % n) as usize
    }

    fn cell_id(ix: usize, iy: usize, iz: usize, n: usize) -> usize {
        (ix * n + iy) * n + iz
    }

    //fn cell_id_to_3d

    fn position_to_cell_3d(
        pos: &Vector3<f64>,
        box_size: &Vector3<f64>,
        n_cells: usize,
    ) -> (usize, usize, usize) {
        /*
        Map to [0, L) first if you're storing unwrapped coordinates
        if your pos is already in [0, L), you can skip this
         */

        let x = pos.x.rem_euclid(box_size.x); // wrap around x coordinate
        let y = pos.y.rem_euclid(box_size.y); // wrap around y coordinate
        let z = pos.z.rem_euclid(box_size.z); // wrap around z coordinate

        // normalized coordinates around [0, 1)

        let fx = (x / box_size.x) * n_cells as f64;
        let fy = (y / box_size.y) * n_cells as f64;
        let fz = (z / box_size.z) * n_cells as f64;

        // floor gives 0..n_cells-1, but we still wrap to be safe
        let ix = wrap_index(fx.floor() as isize, n_cells);
        let iy = wrap_index(fy.floor() as isize, n_cells);
        let iz = wrap_index(fz.floor() as isize, n_cells);

        (ix, iy, iz)
    }

    pub struct SimulationBox {
        pub x_dimension: f64,
        pub y_dimension: f64,
        pub z_dimension: f64,
    }

    pub struct MolecularCoordinates {
        /*
        struct to store information about each subdivsion of cells
        */
        pub center: Vector3<f64>,
        //pub length: f64,
        pub half_length: Vector3<f64>,
        pub index: Vector3<usize>,
        pub atom_index: Vec<usize>,
        pub head: Vec<Option<usize>>,
        pub next: Vec<Option<usize>>,
    }

    impl SimulationBox {
        pub fn create_subcells(&self, n_cells: usize) -> Vec<MolecularCoordinates> {
            let n_cells_per_dim = n_cells;
            let n_cells_total = n_cells_per_dim * n_cells_per_dim * n_cells_per_dim;
            let mut cells = Vec::with_capacity(n_cells_total); // create the cells

            let dx = self.x_dimension / n_cells_per_dim as f64;
            let dy = self.y_dimension / n_cells_per_dim as f64;
            let dz = self.z_dimension / n_cells_per_dim as f64;

            for ix in 0..n_cells_per_dim {
                for iy in 0..n_cells_per_dim {
                    for iz in 0..n_cells_per_dim {
                        // create the cells and push the coordinates (with the center implemented)
                        cells.push(MolecularCoordinates {
                            center: Vector3::new(
                                (ix as f64 + 0.5) * dx,
                                (iy as f64 + 0.5) * dy,
                                (iz as f64 + 0.5) * dz,
                            ),
                            half_length: Vector3::new(dx * 0.5, dy * 0.5, dz * 0.5),
                            index: Vector3::new(ix, iy, iz),
                            atom_index: Vec::new(),
                            head: vec![None; n_cells_total],
                            next: Vec::new(), //
                        });
                    }
                }
            }

            cells
        }

        pub fn store_atoms_in_cells_particles(
            &self,
            particles: &mut Vec<Particle>,
            cells: &mut Vec<MolecularCoordinates>, // the created cells
            n_cells: usize,
        ) -> () {
            /*
            I don't need to store all the coordinates, just the indices of the atoms/systems
             */

            // Clear previous contents
            for cell in cells.iter_mut() {
                cell.atom_index.clear();
            }

            let box_size = Vector3::new(self.x_dimension, self.y_dimension, self.z_dimension);
            // let n_cells_total = cells.len();

            // Store all the particles in the appropriate

            for (i, particle) in particles.iter_mut().enumerate() {
                let (ix, iy, iz) = position_to_cell_3d(&particle.position, &box_size, n_cells);

                let cid = cell_id(ix, iy, iz, n_cells);
                cells[cid].atom_index.push(i);
            }
        }

        pub fn store_atoms_in_cells_systems(
            &self,
            systems: &mut Vec<System>,
            cells: &mut Vec<MolecularCoordinates>,
            n_cells: usize,
        ) -> () {
            /*
            I don't need to store all the coordinates, just the indices of the atoms/systems
             */

            // Clear previous contents
            for cell in cells.iter_mut() {
                cell.atom_index.clear();
            }

            let box_size = Vector3::new(self.x_dimension, self.y_dimension, self.z_dimension);

            for (i, system) in systems.iter_mut().enumerate() {
                for particle in system.atoms.iter() {
                    let (ix, iy, iz) = position_to_cell_3d(&particle.position, &box_size, n_cells);
                    let cid = cell_id(ix, iy, iz, n_cells);
                    cells[cid].atom_index.push(i);
                }
            }
        }
    }
}

pub mod lennard_jones_simulations {

    use super::*; //

    use crate::cell::cell::{CellList, Vec3};
    use crate::error::error::compute_average_val;
    use crate::parameters::lj_parameters::lennard_jones_potential;
    use crate::thermostat_barostat::andersen::andersen::apply_andersen_collisions;
    use crate::thermostat_barostat::nose_hoover::nose_hoover::apply_thermostat_nose_hoover_particles;

    use log::{debug, info};
    #[cfg(feature = "mpi")]
    use mpi::collective::SystemOperation;
    #[cfg(feature = "mpi")]
    use mpi::traits::*;
    use nalgebra::{zero, Vector3};
    use rand::Rng;
    use rand_distr::{Distribution, Normal};

    // importing bonds
    use crate::molecule::molecule::make_h2_system;
    use crate::molecule::molecule::Bond;
    use crate::molecule::molecule::System;
    use crate::molecule::molecule::{
        apply_all_bonded_forces_and_energy, apply_bonded_forces_and_energy, Angle, Dihedral,
        Improper,
    };

    use crate::lennard_jones_simulations::cell_subdivision::MolecularCoordinates;

    #[derive(Clone, Debug)]
    pub struct LJParameters {
        // lennard jones parameters and the number of atoms that we have of that parameter
        pub epsilon: f64,
        pub sigma: f64,
        pub number_of_atoms: i32,
    }

    #[derive(Clone, Debug)]
    pub struct Particle {
        pub id: usize,
        pub position: Vector3<f64>,
        pub velocity: Vector3<f64>,
        pub force: Vector3<f64>,
        pub lj_parameters: LJParameters,
        pub mass: f64,
        pub energy: f64,
        pub atom_type: f64,
        pub charge: f64,
    }

    #[derive(Clone)]
    pub struct SimulationSummary {
        pub energy: f64,
    }

    pub enum InitOutput {
        Particles(Vec<Particle>), // define a particles system (single point particle)
        Systems(Vec<System>),     // Define actual molecules
    }

    pub enum InitMode {
        Atoms,
        Molecules,
    }

    impl Particle {
        pub fn distance(&self, other: &Particle) -> f64 {
            // Compute the distance between two particles
            (self.position - other.position).norm()
        }

        pub fn maxwellboltzmannvelocity(&mut self, temp: f64, mass: f64, _v_max: f64) {
            let mut rng = rand::rng();

            // Standard deviation based on MB distribution
            let sigma_mb = (temp / mass).sqrt();

            // Create a normal distribution with mean = 0, std = sigma
            let normal = Normal::new(0.0, sigma_mb).unwrap();

            // Assign each velocity component independently
            self.velocity[0] = normal.sample(&mut rng);
            self.velocity[1] = normal.sample(&mut rng);
            self.velocity[2] = normal.sample(&mut rng);

            debug!("Assigned Maxwell-Boltzmann velocity: {:?}", self.velocity);
        }

        fn update_position_verlet(&mut self, dt: f64) -> () {
            /*
            Verlet scheme to change the position
            Use the verlet scheme to change the velocity
            */
            //let a = self.force / self.mass;
            self.position += self.velocity * dt; //+ 0.5 * a * dt * dt;
        }

        fn update_velocity_verlet(&mut self, accelerations: Vector3<f64>, dt: f64) {
            /*
            Verlet scheme to update the velocity
             */
            self.velocity += 0.5 * accelerations * dt;
        }
    }

    pub fn site_site_energy_calculation(particles: &mut Vec<Particle>, box_length: f64) -> f64 {
        /*
        Computing the total Lennard-Jones energy between all distinct pairs of particles in a molecular system,
        using site-site interactions

        The coordinates r_ia of a site a in molecule i are stored in the elements r(:, i, a)

            For example, if we have two diatomic molecules, then we have r_1a (site a of molecule 1) and
            r_2a (site a of molecule 2). Each molecule is a diatomic molecule (for example, O2).

            We already have a set of particles with the lennard jones parameters defined and stored within. Using
            that data, we need to compute the site_site energy

         */

        let mut total_energy = 0.0;

        for i in 0..particles.len() {
            for j in (i + 1)..particles.len() {
                // double loop over all coordinates in the system

                let sigma_i = particles[i].lj_parameters.sigma; // for particle i, get the sigma
                let epsilon_i = particles[i].lj_parameters.epsilon; // for particle i, get the epsilon
                let sigma_j = particles[j].lj_parameters.sigma; // for particle j, get the sigma
                let epsilon_j = particles[j].lj_parameters.epsilon; // for particle j, get the epsilon

                // Using Lorentz-Bethelot mixing rules
                let computed_sigma = (sigma_i + sigma_j) / 2.0;
                let computed_epsilon = (epsilon_i * epsilon_j).sqrt();
                let r_vec = particles[j].position - particles[i].position; // We have already applied PBC to wrap the positions
                let r_vec_mic = minimum_image_convention(r_vec, box_length); // minimum image convention is used for computing the true closest distance between the partcle i and j, through the images rather than take the longest distance from within the same image
                let r = r_vec_mic.norm();
                let potential = lennard_jones_potential(r, computed_sigma, computed_epsilon);

                // Sum the total energy with the pairwise potential in the system
                total_energy += potential;
            }
        }
        total_energy
    }

    pub fn set_molecular_positions_and_velocities(
        system_mol: &mut InitOutput,
        temp: f64,
        box_length: f64,
    ) -> () {
        let mut rng = rand::rng();
        // loop over each moleucle

        // Create a normal distribution with mean = 0, std = sigma

        match system_mol {
            InitOutput::Particles(particles) => {
                for index in 0..particles.len() {
                    // Each element is a System
                    particles[index].velocity[0] = rng.random_range(-1.0..1.0);
                    particles[index].velocity[1] = rng.random_range(-1.0..1.0);
                    particles[index].velocity[2] = rng.random_range(-1.0..1.0);
                }
            }
            InitOutput::Systems(systems) => {
                for sys in systems.iter_mut() {
                    // Randomize each molecule position as a rigid translation so cloned
                    // systems do not remain perfectly overlapped at initialization.
                    let translation = Vector3::new(
                        rng.random_range(0.0..box_length),
                        rng.random_range(0.0..box_length),
                        rng.random_range(0.0..box_length),
                    );

                    // Each element is a System
                    // loop over each atom
                    for atom in sys.atoms.iter_mut() {
                        atom.position += translation;
                        let sigma_mb = (temp / atom.mass).sqrt(); // 1.0 needs to be replaced with mass
                        let normal = Normal::new(0.0, sigma_mb).unwrap();

                        atom.velocity[0] = normal.sample(&mut rng);
                        atom.velocity[1] = normal.sample(&mut rng);
                        atom.velocity[2] = normal.sample(&mut rng);
                    }
                }
            }
        }
    }

    pub fn create_atoms_with_set_positions_and_velocities(
        number_of_atoms: i64,
        temp: f64,
        mass: f64,
        v_max: f64,
        box_dim_max: f64,
        use_atom: bool,
    ) -> Result<InitOutput, String> {
        /*

        Create N atoms with temperature and mass

        Here, we are going with the assumption that we are creating a simulation box
        that is cubic, meaning that we will get a (min + max) * (min + max) *  (min+max) volume
        system for all the molecules

         */
        if box_dim_max <= 0.0 {
            return Err("box_dim_max must be > 0.0".to_string());
        }

        let mut vector_positions: Vec<Particle> = Vec::new();
        let mut vector_system_positions: Vec<System> = Vec::new();
        let mut rng = rand::rng();
        // Create the number of atoms in the system with the system as necessary
        if !use_atom {
            for index in 0..number_of_atoms {
                // For each particle, we wish to create the initial posiiton,
                // the velocity, and the LJ parameters attached to it
                let mut particle = Particle {
                    // create position for the atom in question
                    position: Vector3::new(
                        // generate x y z position values using the configured box size
                        rng.random_range(0.0..box_dim_max),
                        rng.random_range(0.0..box_dim_max),
                        rng.random_range(0.0..box_dim_max),
                    ),

                    // create velocity for atom in question
                    velocity: Vector3::new(
                        // generate velocity values between -1 and 1
                        rng.random_range(-1.0..1.0),
                        rng.random_range(-1.0..1.0),
                        rng.random_range(-1.0..1.0),
                    ),
                    // Create the LJ parameters for the atom - default parameters are 1.0, 1.0
                    lj_parameters: (LJParameters {
                        epsilon: 1.0,
                        sigma: 1.0,
                        number_of_atoms: 3,
                    }),
                    force: zero(), // initial force on the atom
                    mass: mass,    // the mass
                    energy: 0.0,
                    atom_type: 0.0,
                    charge: 0.0,
                    id: index as usize,
                };

                // Reset the positions to the maxwell boltzmann distibution of velocities
                particle.maxwellboltzmannvelocity(temp, mass, v_max);
                // push those values into the vector
                vector_positions.push(particle); // push the newly assigned particle into the positions
            }
            Ok(InitOutput::Particles(vector_positions))
        } else {
            // This needs to be fixed
            for _ in 0..number_of_atoms {
                let h2_system = make_h2_system(); //
                vector_system_positions.push(h2_system);
            }
            Ok(InitOutput::Systems(vector_system_positions))
        }
    }

    pub fn implement_shake() -> () {}

    //pub fn run_verlet_update_nve(state: &mut InitOutput, dt: f64, box_length: f64) -> () {
    //    /*
    //    Update the position and velocity of the particle using the verlet scheme
    //     */
    //    match state {
    //        InitOutput::Particles(particles) => {
    //            for particle in particles.iter_mut() {
    //                particle.update_position_verlet(dt);
    //            }
    //
    //            pbc_update(particles, box_length);
    //
    //            let simulation_box = cell_subdivision::SimulationBox {
    //                x_dimension: box_length,
    //                y_dimension: box_length,
    //                z_dimension: box_length,
    //            };
    //
    //            // Create the subcells - here we have used a subdivision of 10 for the cells
    //            let mut subcells = simulation_box.create_subcells(10);
    //
    //            // Store the coordinates in cells
    //            simulation_box.store_atoms_in_cells_particles(particles, &mut subcells, 10);
    //            compute_forces_particles(particles, box_length, &mut subcells);
    //
    //            for particle in particles.iter_mut() {
    //                println!(
    //                    "The original position and velocity is {:?} and {:?} ",
    //                    particle.position, particle.velocity
    //                );
    //                let a_new = particle.force / particle.mass;
    //                particle.update_velocity_verlet(a_new, dt);
    //
    //                println!(
    //                    "After a iteration step, the position and velocity is {:?} and {:?} ",
    //                    particle.position, particle.velocity
    //                );
    //            }
    //        }
    //
    //        // for each 'system' - actual molecule in the simulation
    //        InitOutput::Systems(systems) => {
    //            for sys in systems.iter_mut() {
    //                for s in sys.atoms.iter_mut() {
    //                    // this is a vector of particles
    //                    println!(
    //                        "The original position and velocity is {:?} and {:?} for the system",
    //                        s.position, s.velocity
    //                    );
    //                    s.update_position_verlet(dt);
    //                }
    //                pbc_update(&mut sys.atoms, box_length);
    //            }
    //
    //            let simulation_box = cell_subdivision::SimulationBox {
    //                x_dimension: box_length,
    //                y_dimension: box_length,
    //                z_dimension: box_length,
    //            };
    //
    //            // Create the subcells - here we have used a subdivision of 10 for the cells
    //            let mut subcells = simulation_box.create_subcells(10);
    //            // Store the coordinates in cells
    //            simulation_box.store_atoms_in_cells_systems(systems, &mut subcells, 10);
    //
    //            // update the velocity
    //            for sys in systems.iter_mut() {
    //                for s in sys.atoms.iter_mut() {
    //                    let a_new = s.force / s.mass; // compu
    //                    s.update_velocity_verlet(a_new, dt);
    //
    //                    println!(
    //                        "After a iteration step, the position and velocity is {:?} and {:?} ",
    //                        s.position, s.velocity
    //                    );
    //                }
    //            }
    //        }
    //    }
    //}

    pub fn apply_bond_force(particles: &mut [Particle], b: &Bond, box_length: f64) -> f64 {
        let rij = particles[b.atom1].position - particles[b.atom2].position;
        let rij_mic = minimum_image_convention(rij, box_length);
        let r = safe_norm(rij_mic.norm());
        let dr = r - b.r0;
        let f_mag = -b.k * dr; // along r̂, attractive if r>r0
        let f_vec = (rij_mic / r) * f_mag; // vector force on i

        particles[b.atom1].force += f_vec;
        particles[b.atom2].force -= f_vec;

        0.5 * b.k * dr * dr
    }

    pub fn compute_pair_forces_vector(dr: Vec3, r2: f64, sigma: f64, epsilon: f64) -> Vec3 {
        if r2 == 0.0 {
            return Vec3::zero();
        }
        let r = r2.sqrt();
        // make sure the mixing rules are correct
        let f_mag = lennard_jones_force_scalar(r, sigma, epsilon);
        let f_vec = (dr / r) * f_mag; // along r-hat
        f_vec
    }

    pub fn compute_forces_particles(
        particles: &mut Vec<Particle>,
        box_length: f64,
        cells: &mut Vec<MolecularCoordinates>,
    ) {
        /*
                Computing forces between the single point particles

            Todo - need to add a component where I can only need to loop through
            the cells that are closest and only compute the forces between
        the particles in these cells

                 */
        for p in particles.iter_mut() {
            p.force = Vector3::zeros();
        }
        let mut cell_match_storage: Vec<Vec<i32>> = Vec::new();

        // initalize zero forces for each particle

        // Instead of looping over each particle,
        // will probably need to loop over the
        // cells that match within the distance boundary
        let cell_limit = 1.0;

        for i in 0..cells.len() {
            for j in (i + 1)..cells.len() {
                // compute the distance between the cells
                let cell_dist = cells[i].center - cells[j].center;
                let cell_r = cell_dist.norm();
                if cell_r <= cell_limit {
                    let cell_match_instance = vec![i as i32, j as i32];
                    //println!(
                    //    " We have the following cell match instance {:?}",
                    //    cell_match_instance
                    //);
                    cell_match_storage.push(cell_match_instance);
                }
            }
        }
        // dedup the indices
        dedup_permutation(&mut cell_match_storage); // TODO - understand the operation of this deduping

        // loop over the cells
        for cell_i in 0..cell_match_storage.len() {
            // select cell entry
            let cell_i_idx = cell_match_storage[cell_i][0];
            let cell_ii_idx = cell_match_storage[cell_i][1];

            for i_index in cells[cell_i_idx as usize].atom_index.iter() {
                for j_index in cells[cell_ii_idx as usize].atom_index.iter() {
                    let r_vec = particles[*i_index].position - particles[*j_index].position;
                    let r_mic = minimum_image_convention(r_vec, box_length);
                    let r = r_mic.norm(); // compute the distance
                    if r == 0.0 {
                        continue;
                    }
                    let si = particles[*i_index].lj_parameters.sigma;
                    let ei = particles[*i_index].lj_parameters.epsilon;
                    let sj = particles[*j_index].lj_parameters.sigma;
                    let ej = particles[*j_index].lj_parameters.epsilon;
                    let sigma = 0.5 * (si + sj);
                    let epsilon = (ei * ej).sqrt();
                    let f_mag = lennard_jones_force_scalar(r, sigma, epsilon);
                    let f_vec = (r_mic / r) * f_mag; // along r-hat

                    // action = -reaction
                    particles[*i_index].force -= f_vec;
                    particles[*j_index].force += f_vec;
                    debug!(
                        "Pair force update: i={:?}, j={:?}",
                        particles[*i_index].force, particles[*j_index].force
                    );
                }
            }
        }
    }

    fn compute_bonded_forces(
        atoms: &mut Vec<Particle>,
        bonds: &[Bond],
        angles: &[Angle],
        dihedrals: &[Dihedral],
        impropers: &[Improper],
        box_length: f64,
    ) -> f64 {
        apply_all_bonded_forces_and_energy(atoms, bonds, angles, dihedrals, impropers, box_length)
    }

    fn compute_bonded_energy(
        atoms: &mut Vec<Particle>,
        bonds: &[Bond],
        angles: &[Angle],
        dihedrals: &[Dihedral],
        impropers: &[Improper],
        box_length: f64,
    ) -> f64 {
        // NOTE:
        // `apply_all_bonded_forces_and_energy` updates forces as a side effect.
        // During MD we call this helper for diagnostics only (energy reporting),
        // so we must avoid mutating the active force buffers, otherwise bonded
        // forces are effectively applied twice and the integrator can blow up.
        let mut atoms_for_energy = atoms.clone();
        apply_all_bonded_forces_and_energy(
            &mut atoms_for_energy,
            bonds,
            angles,
            dihedrals,
            impropers,
            box_length,
        )
    }

    pub fn compute_intermolecular_forces_systems(systems: &mut [System], box_length: f64) -> f64 {
        /*
        Compute Lennard-Jones interactions between atoms belonging to different systems.
        Intra-molecular interactions are omitted here and handled by bonded terms.
         */
        let mut total_energy = 0.0;

        for i in 0..systems.len() {
            for j in (i + 1)..systems.len() {
                let (left, right) = systems.split_at_mut(j);
                let sys_i = &mut left[i];
                let sys_j = &mut right[0];

                for atom_i in sys_i.atoms.iter_mut() {
                    for atom_j in sys_j.atoms.iter_mut() {
                        let r_vec = atom_j.position - atom_i.position;
                        let r_mic = minimum_image_convention(r_vec, box_length);
                        let r = safe_norm(r_mic.norm());

                        let sigma = 0.5 * (atom_i.lj_parameters.sigma + atom_j.lj_parameters.sigma);
                        let epsilon =
                            (atom_i.lj_parameters.epsilon * atom_j.lj_parameters.epsilon).sqrt();

                        let f_mag = lennard_jones_force_scalar(r, sigma, epsilon);
                        let f_vec = (r_mic / r) * f_mag;

                        atom_i.force -= f_vec;
                        atom_j.force += f_vec;

                        total_energy += lennard_jones_potential(r, sigma, epsilon);
                    }
                }
            }
        }

        total_energy
    }

    pub fn intermolecular_site_site_energy_systems(systems: &[System], box_length: f64) -> f64 {
        /*
        Compute Lennard-Jones potential energy between atoms in different systems.
         */
        let mut total_energy = 0.0;

        for i in 0..systems.len() {
            for j in (i + 1)..systems.len() {
                let sys_i = &systems[i];
                let sys_j = &systems[j];

                for atom_i in sys_i.atoms.iter() {
                    for atom_j in sys_j.atoms.iter() {
                        let r_vec = atom_j.position - atom_i.position;
                        let r_mic = minimum_image_convention(r_vec, box_length);
                        let r = safe_norm(r_mic.norm());

                        let sigma = 0.5 * (atom_i.lj_parameters.sigma + atom_j.lj_parameters.sigma);
                        let epsilon =
                            (atom_i.lj_parameters.epsilon * atom_j.lj_parameters.epsilon).sqrt();

                        total_energy += lennard_jones_potential(r, sigma, epsilon);
                    }
                }
            }
        }

        total_energy
    }

    pub fn compute_bonded_forces_system(
        atoms: &mut Vec<Particle>,
        bonds: &[Bond],
        box_length: f64,
    ) -> f64 {
        /*
        Initialize the forces on the systems (molecules) in the simulation box, and apply newton's third law to each system (molecule)
         */
        for a in atoms.iter_mut() {
            a.force = Vector3::zeros();
        }
        // compute the intermolecular forces - TODO - tbis should omitt
        // intramolecular forces

        apply_bonded_forces_and_energy(atoms, bonds, box_length)
    }

    pub fn compute_all_bonded_forces_system(system: &mut System, box_length: f64) -> f64 {
        for atom in system.atoms.iter_mut() {
            atom.force = Vector3::zeros();
        }
        compute_bonded_forces(
            &mut system.atoms,
            &system.bonds,
            &system.angles,
            &system.dihedrals,
            &system.impropers,
            box_length,
        )
    }

    pub fn compute_electrostatic_forces_systems(systems: &mut [System], box_length: f64) -> f64 {
        let pme = PmeConfig::default();
        add_electrostatic_forces_systems(systems, box_length, &pme)
    }

    // -- temperature related computations

    pub fn kinetic_energy_particles(particles: &[Particle]) -> f64 {
        particles
            .iter()
            .map(|p| 0.5 * p.mass * p.velocity.norm_squared())
            .sum()
    }

    pub fn compute_temperature_particles(particles: &[Particle], dof: usize) -> f64 {
        const KB_KJ_PER_MOL_K: f64 = 0.008_314_462_618_153_24;
        if dof == 0 {
            return 0.0;
        }
        2.0 * kinetic_energy_particles(particles) / (KB_KJ_PER_MOL_K * dof as f64)
    }

    pub fn compute_temperature(state: &mut InitOutput, dof: usize) -> f64 {
        /*
                    Sums the kinetic energy over the whole system

                K = total (0.5 * m * v ^ 0.5)

            T = 2K / (K_b (dof))

        // TODO - need to actually implement the boltzmann constant for computing the temperature

         */

        let total_ke = match state {
            InitOutput::Particles(p) => kinetic_energy_particles(p),
            InitOutput::Systems(systems) => {
                if dof == 0 {
                    return 0.0;
                }
                systems
                    .iter()
                    .map(|sys| kinetic_energy_particles(&sys.atoms))
                    .sum()
            }
        };

        const KB_KJ_PER_MOL_K: f64 = 0.008_314_462_618_153_24;
        2.0 * total_ke / (KB_KJ_PER_MOL_K * dof as f64)
    }

    pub fn compute_pressure_particles(particles: &[Particle], box_length: f64) -> f64 {
        let n = particles.len();
        if n == 0 || box_length <= 0.0 {
            return 0.0;
        }

        let dof = 3 * n; // each particle has a degree of freedom of 3
        let temperature = compute_temperature_particles(particles, dof);
        let volume = box_length.powi(3);

        // The virial measures how particle positions correlate with force

        // hence - virial shows how strongly forces act at a given distance

        let mut virial = 0.0;

        for i in 0..n {
            // loop over the particles
            for j in (i + 1)..n {
                let r_vec = particles[j].position - particles[i].position; // compute the distance between particles
                let r_mic = minimum_image_convention(r_vec, box_length); // recompute the distance according to the minimum image convention
                let r = r_mic.norm(); // gets the magnitude of the distance

                if r == 0.0 {
                    continue;
                }

                let si = particles[i].lj_parameters.sigma; // sigma for the lennard jones force contribution
                let ei = particles[i].lj_parameters.epsilon; // epsilon parameter

                let sj = particles[j].lj_parameters.sigma;
                let ej = particles[j].lj_parameters.epsilon;

                // need to remind myself which mixing rules are these
                let sigma = 0.5 * (si + sj);
                let epsilon = (ei * ej).sqrt();

                let f_mag = lennard_jones_force_scalar(r, sigma, epsilon); // pairwise force
                let f_vec = (r_mic / r) * f_mag; //we get the unit vector, and we get the radial force magnitude

                virial += r_mic.dot(&f_vec);
            }
        }
        /*
        Compute instantaneous pressure using the virial equation of state.

        P = (N T)/V + (1/(3V)) * Σ_{i<j} (r_ij · F_ij)

        where:
        - N      = number of particles
        - T      = instantaneous temperature (from kinetic energy via equipartition)
        - V      = simulation box volume
        - r_ij   = minimum-image displacement vector between particles i and j
        - F_ij   = force on particle i due to particle j
        - Σ_{i<j} ensures each pair interaction is counted once

        Interpretation:
        - First term: ideal gas contribution from kinetic motion
        P_ideal = N T / V

        - Second term: configurational (virial) contribution from interactions
        P_virial = (1/(3V)) Σ r_ij · F_ij

        Assumptions:
        - Reduced units with k_B = 1
        - Pairwise-additive forces (Lennard-Jones here)
        - Periodic boundary conditions using minimum image convention

        Total pressure:
         */
        (n as f64 * temperature) / volume + virial / (3.0 * volume)
    }

    pub fn apply_thermostat(state: &mut InitOutput, target_temperature: f64) {
        match state {
            InitOutput::Particles(particles) => {
                // dof: subtract 3 to account for removing COM motion (classic MD trick)
                let dof = 3 * particles.len();
                if dof == 0 {
                    return;
                }

                let current_temperature = compute_temperature_particles(particles, dof);
                if current_temperature == 0.0 {
                    return;
                }

                let lambda = (target_temperature / current_temperature).sqrt();

                for p in particles.iter_mut() {
                    p.velocity *= lambda;
                }

                info!(
                    "Thermostat[particles] T: {:.2} -> {:.2} (scale λ={:.4})",
                    current_temperature, target_temperature, lambda
                );
            }

            InitOutput::Systems(systems) => {
                // Option A: rescale each system independently (simple & clear)
                for (si, sys) in systems.iter_mut().enumerate() {
                    let natoms = sys.atoms.len();
                    if natoms == 0 {
                        continue;
                    }

                    let dof = 3 * natoms;
                    if dof == 0 {
                        continue;
                    }

                    let current_temperature = compute_temperature_particles(&sys.atoms, dof);
                    if current_temperature == 0.0 {
                        continue;
                    }

                    let lambda = (target_temperature / current_temperature).sqrt();

                    for a in sys.atoms.iter_mut() {
                        a.velocity *= lambda;
                    }

                    info!(
                        "Thermostat[system {si}] T: {:.2} -> {:.2} (scale λ={:.4})",
                        current_temperature, target_temperature, lambda
                    );
                }

                // Option B (alternative): compute one global T over all atoms, apply single λ.
                // Implement later if you want physically consistent global NVT.
            }
        }
    }

    pub fn apply_thermostat_berendsen_particles(
        particles: &mut Vec<Particle>,
        target_temperature: f64,
        tau: f64,
        dt: f64,
    ) -> () {
        /*

        If the system's instantanepus temperature T differs from the target temperature T_0, the Berendsen
        thermostat weakly couples the system to a head bath that gently nudges T towards T_0 over a characteristic
        relaxation time

         */
        let dof = 3 * particles.len();
        let current_temperature = compute_temperature_particles(particles, dof);
        // Bail if parameters are nonsense or temperature is zero/negative
        if tau <= 0.0 || dt <= 0.0 || current_temperature <= 0.0 || target_temperature <= 0.0 {
            return;
        }

        // Discrete Berendsen: T' = T * (1 + (dt/tau)(T_0/T - 1))
        // Velocitiea s scale as sqrt (T'/T)
        let x = (dt / tau) * (target_temperature / current_temperature - 1.0);

        // clamp to avoid negative
        let x_clamped = x.clamp(-0.9, 10.0);
        let lambda = (1.0 + x_clamped).max(1e-12).sqrt();

        for particle in particles {
            particle.velocity *= lambda;
        }
    }

    pub fn apply_thermostat_andersen_particles(
        particles: &mut Vec<Particle>,
        box_length: f64,
        target_temperature: f64,
        dt: f64,
        t_max: f64,
    ) -> () {
        /*
        Initialize system and compute the forces and energy
         */
        let mut t = 0.0;
        let mut switch = 1;

        while t < t_max {
            // Propagates the half step
            run_md_andersen_particles(particles, dt, box_length, target_temperature, 1.0, switch);

            let simulation_box = cell_subdivision::SimulationBox {
                x_dimension: box_length,
                y_dimension: box_length,
                z_dimension: box_length,
            };

            // Create the subcells - here we have used a subdivision of 10 for the cells
            let mut subcells = simulation_box.create_subcells(10);
            // Store the coordinates in cells
            simulation_box.store_atoms_in_cells_particles(particles, &mut subcells, 10);
            // Compute the forces in the system
            compute_forces_particles(particles, box_length, &mut subcells);
            // switches to 2
            switch = 2;
            // Propagates the second half time step
            run_md_andersen_particles(particles, dt, box_length, target_temperature, 1.0, switch);
            t = t + dt;
        }
    }

    pub fn apply_thermostat_berendsen(
        state: &mut InitOutput,
        target_temperature: f64,
        tau: f64,
        dt: f64,
    ) {
        match state {
            InitOutput::Particles(particles) => {
                apply_thermostat_berendsen_particles(particles, target_temperature, tau, dt);
            }
            InitOutput::Systems(systems) => {
                // Option A: per-molecule coupling
                for sys in systems.iter_mut() {
                    apply_thermostat_berendsen_particles(
                        &mut sys.atoms,
                        target_temperature,
                        tau,
                        dt,
                    );
                }

                // Option B (if you prefer one global T and λ across all atoms):
                //  - flatten all atoms, compute global T, single λ
                //  - apply to every atom in every sys
                // Do that later if you care about strict ensemble correctness.
            }
        }
    }

    pub fn pbc_update(particles: &mut Vec<Particle>, box_length: f64) {
        /*
        Depending on what kind of system we are injecting to this function, we want to produce the correct
        pbc update to the coordinates
         */
        for particle in particles.iter_mut() {
            for i in 0..3 {
                particle.position[i] = particle.position[i].rem_euclid(box_length);
            }
        }
    }

    pub fn compute_total_energy_and_print(state: &mut InitOutput, box_length: f64) -> f64 {
        /*
        compute the total kinetic + potential energy of the system
         */
        let mut kinetic_energy = 0.0;
        let mut potential_energy = 0.0;

        match state {
            InitOutput::Particles(particles) => {
                for p in particles.iter_mut() {
                    let v2 = p.velocity.norm_squared();
                    kinetic_energy += 0.5 * p.mass * v2;
                }
                potential_energy = site_site_energy_calculation(particles, box_length);
            }

            InitOutput::Systems(systems) => {
                for sys in systems.iter_mut() {
                    for a in sys.atoms.iter() {
                        let v2 = a.velocity.norm_squared();
                        kinetic_energy += 0.5 * a.mass * v2;
                    }
                    potential_energy = site_site_energy_calculation(&mut sys.atoms, box_length);
                }
            }
        }

        debug!(
            "Total instantaneous energy: {:.6}",
            kinetic_energy + potential_energy
        );

        kinetic_energy + potential_energy
    }

    pub fn minimum_image_convention(rij: Vector3<f64>, box_length: f64) -> Vector3<f64> {
        Vector3::new(
            rij[0] - box_length * (rij[0] / box_length).round(),
            rij[1] - box_length * (rij[1] / box_length).round(),
            rij[2] - box_length * (rij[2] / box_length).round(),
        )
    }

    fn erfc_approx(x: f64) -> f64 {
        // Abramowitz and Stegun 7.1.26
        let z = x.abs();
        let t = 1.0 / (1.0 + 0.3275911 * z);
        let a1 = 0.254829592;
        let a2 = -0.284496736;
        let a3 = 1.421413741;
        let a4 = -1.453152027;
        let a5 = 1.061405429;
        let poly = (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t;
        let erf = 1.0 - poly * (-z * z).exp();
        let erf = if x < 0.0 { -erf } else { erf };
        1.0 - erf
    }

    fn add_electrostatic_real_space_particles(
        particles: &mut [Particle],
        box_length: f64,
        pme: &PmeConfig,
    ) -> f64 {
        let mut energy = 0.0;
        let alpha = pme.alpha;
        let rc = pme.real_cutoff;
        let k_e = coulomb_prefactor();

        for i in 0..particles.len() {
            for j in (i + 1)..particles.len() {
                let qi = particles[i].charge;
                let qj = particles[j].charge;
                if qi == 0.0 && qj == 0.0 {
                    continue;
                }

                let rij = minimum_image_convention(
                    particles[j].position - particles[i].position,
                    box_length,
                );
                let r = rij.norm();
                if r <= 1e-12 || r > rc {
                    continue;
                }

                let ar = alpha * r;
                let erfc_ar = erfc_approx(ar);
                let exp_term = (-(ar * ar)).exp();
                let qq = k_e * qi * qj;

                energy += qq * erfc_ar / r;

                let scalar = qq
                    * (erfc_ar / (r * r)
                        + (2.0 * alpha / std::f64::consts::PI.sqrt()) * exp_term / r);
                let f_vec = -(rij / r) * scalar;
                particles[i].force += f_vec;
                particles[j].force -= f_vec;
            }
        }

        energy
    }

    fn add_electrostatic_reciprocal_particles(
        particles: &mut [Particle],
        box_length: f64,
        pme: &PmeConfig,
    ) -> f64 {
        if particles.is_empty() {
            return 0.0;
        }

        let volume = box_length.powi(3);
        let alpha = pme.alpha;
        let kmax = pme.kmax;
        let k_e = coulomb_prefactor();
        let two_pi_over_l = 2.0 * std::f64::consts::PI / box_length;
        let mut energy = 0.0;

        for nx in -kmax..=kmax {
            for ny in -kmax..=kmax {
                for nz in -kmax..=kmax {
                    if nx == 0 && ny == 0 && nz == 0 {
                        continue;
                    }

                    let kvec = Vector3::new(
                        nx as f64 * two_pi_over_l,
                        ny as f64 * two_pi_over_l,
                        nz as f64 * two_pi_over_l,
                    );
                    let k2 = kvec.norm_squared();
                    if k2 <= 1e-12 {
                        continue;
                    }

                    let damp = (-k2 / (4.0 * alpha * alpha)).exp();
                    let coef = (2.0 * std::f64::consts::PI * k_e / volume) * damp / k2;

                    let mut s_cos = 0.0;
                    let mut s_sin = 0.0;
                    for p in particles.iter() {
                        let phase = kvec.dot(&p.position);
                        s_cos += p.charge * phase.cos();
                        s_sin += p.charge * phase.sin();
                    }

                    energy += coef * (s_cos * s_cos + s_sin * s_sin);

                    for p in particles.iter_mut() {
                        let phase = kvec.dot(&p.position);
                        let sin_i = phase.sin();
                        let cos_i = phase.cos();
                        let force_coeff =
                            -(4.0 * std::f64::consts::PI * k_e * p.charge / volume) * damp / k2;
                        let proj = s_cos * sin_i - s_sin * cos_i;
                        p.force += kvec * (force_coeff * proj);
                    }
                }
            }
        }

        let self_energy: f64 = particles.iter().map(|p| p.charge * p.charge).sum::<f64>()
            * (-k_e * alpha / std::f64::consts::PI.sqrt());

        energy + self_energy
    }

    fn add_electrostatic_forces_particles(
        particles: &mut [Particle],
        box_length: f64,
        pme: &PmeConfig,
    ) -> f64 {
        add_electrostatic_real_space_particles(particles, box_length, pme)
            + add_electrostatic_reciprocal_particles(particles, box_length, pme)
    }

    fn add_electrostatic_forces_systems(
        systems: &mut [System],
        box_length: f64,
        pme: &PmeConfig,
    ) -> f64 {
        let mut all_atoms: Vec<Particle> = systems
            .iter()
            .flat_map(|s| s.atoms.iter().cloned())
            .collect();

        let energy = add_electrostatic_forces_particles(&mut all_atoms, box_length, pme);

        let mut idx = 0usize;
        for sys in systems.iter_mut() {
            for atom in sys.atoms.iter_mut() {
                atom.force += all_atoms[idx].force;
                idx += 1;
            }
        }

        energy
    }

    pub fn run_md_andersen_particles(
        particles: &mut Vec<Particle>,
        dt: f64,
        box_length: f64,
        temp: f64,
        nu: f64, // this is the collision frequency
        switch: i64,
    ) -> () {
        // Equations of motion - Andersen thermostat
        let mut a_old: Vec<Vector3<f64>> = Vec::with_capacity(particles.len());
        for a in particles.iter() {
            a_old.push(a.force / a.mass); // compute the acceleration
        }

        if switch == 1 {
            for (p, a_o) in particles.iter_mut().zip(a_old.iter()) {
                // first step velocity verlet
                p.position += dt * p.velocity + (dt * dt) * a_o / 2.0; // update the position current time
                p.velocity += 0.5 * a_o * dt; // update velocity
            }
        } else if switch == 2 {
            /*
            Forces should be recomputed BEFORE this half-kikc
             */

            let simulation_box = cell_subdivision::SimulationBox {
                x_dimension: box_length,
                y_dimension: box_length,
                z_dimension: box_length,
            };

            let mut subcells = simulation_box.create_subcells(10);
            // Store the coordinates in cells
            simulation_box.store_atoms_in_cells_particles(particles, &mut subcells, 10);
            compute_forces_particles(particles, box_length, &mut subcells);

            for p in particles.iter_mut() {
                let a_new = p.force / p.mass; // compute the new acceleration
                p.velocity += 0.5 * a_new * dt;
            }

            /*
            Andersen thermostat step - randomize some velocities.
            probability ~ nu * dt ( valid when nu * dt is small; otherwise use 1 - exp(-nu * dt)
             */

            let mut rng = rand::rng();
            let p_coll = nu * dt;

            for p in particles.iter_mut() {
                let r: f64 = rng.random(); // randomly assign a value between 0 and 1
                if r < p_coll {
                    // If the value of r is smaller than p_col, we reassign the velocity according
                    // to the maxwell boltzmann distribution of that temperature
                    let sigma = (temp / p.mass).sqrt();
                    let normal = Normal::new(0.0, sigma).expect("blah");

                    // assign the new velocities
                    p.velocity = Vector3::new(
                        normal.sample(&mut rng),
                        normal.sample(&mut rng),
                        normal.sample(&mut rng),
                    );
                }
            }
        }
    }

    pub fn run_md_nve_particles(
        particles: &mut Vec<Particle>,
        number_of_steps: i32,
        dt: f64,
        box_length: f64,
        thermostat: &str,
        cutoff: f64,
    ) {
        let mut values: Vec<f32> = Vec::new();
        let box_len = Vec3::new(box_length, box_length, box_length);
        let mut cl = CellList::new(box_len, cutoff); // TODO CELL
        let pme = PmeConfig::default();

        // Create the subcells - we need to have a initial force list for each cell
        cl.rebuild(&particles);
        let mut initial_forces = vec![Vector3::<f64>::zeros(); particles.len()];
        cl.for_each_neighbor_pair(particles, |i, j, dr, r2| {
            let si = particles[i].lj_parameters.sigma;
            let ei = particles[i].lj_parameters.epsilon;
            let sj = particles[j].lj_parameters.sigma;
            let ej = particles[j].lj_parameters.epsilon;
            let sigma = 0.5 * (si + sj);
            let epsilon = (ei * ej).sqrt();

            let f_vec = compute_pair_forces_vector(dr, r2, sigma, epsilon);
            let fv = Vector3::new(f_vec.x, f_vec.y, f_vec.z);

            initial_forces[i] -= fv;
            initial_forces[j] += fv;
        });

        for (p, f) in particles.iter_mut().zip(initial_forces.into_iter()) {
            p.force = f;
        }
        let electrostatic_init_energy =
            add_electrostatic_forces_particles(particles, box_length, &pme);

        let mut kinetic_energy = 0.0;

        // accumulate kinetic energy
        for p in particles.iter() {
            kinetic_energy += 0.5 * p.mass * p.velocity.norm_squared();
        }

        let mut potential_energy =
            site_site_energy_calculation(particles, box_length) + electrostatic_init_energy;
        let mut total_energy = kinetic_energy + potential_energy;

        info!(
            "Init particle energy | E_kin={kinetic_energy:.6} E_pot={potential_energy:.6} E_tot={total_energy:.6}"
        );

        // only used if we are using nose_hoover
        let mut xi_nose_hoover = 0.0;

        // --- time integration loop ---
        for _step in 0..number_of_steps {
            // 1) position updatelet mut forces = vec![Vector3::zeros(); particles.len()];
            let mut forces = vec![Vector3::<f64>::zeros(); particles.len()];

            let mut a_old: Vec<Vector3<f64>> = Vec::with_capacity(particles.len());

            // current acceleration
            for a in particles.iter() {
                a_old.push(a.force / a.mass);
            }

            // 1) velocity update (Verlet - half step)
            for (atom, a_o) in particles.iter_mut().zip(a_old.iter()) {
                atom.velocity += 0.5 * a_o * dt;
            }

            // 2) position update - needs to use the velocity that has been
            // updated using the half step method
            for atom in particles.iter_mut() {
                atom.update_position_verlet(dt);
            }

            // 3) PBC
            pbc_update(particles, box_length);

            cl.rebuild(&particles);

            let mut count = 0usize;

            cl.for_each_neighbor_pair(&particles, |i, j, dr, r2| {
                count += 1;
                // Insert your force calc here (LJ etc.)
                //println!(
                //    "pair ({i},{j}) r2={r2:.4} dr=({:.3},{:.3},{:.3})",
                //    dr.x, dr.y, dr.z
                //);
                //compute_forces_particles_index(&particles, box_length, i as u64, j as u64);

                let si = particles[i as usize].lj_parameters.sigma;
                let ei = particles[i as usize].lj_parameters.epsilon;
                let sj = particles[j as usize].lj_parameters.sigma;
                let ej = particles[j as usize].lj_parameters.epsilon;
                let sigma = 0.5 * (si + sj);
                let epsilon = (ei * ej).sqrt();

                let f_vec = compute_pair_forces_vector(dr, r2, sigma, epsilon);
                let fv = Vector3::new(f_vec.x, f_vec.y, f_vec.z);

                forces[i] -= fv;
                forces[j] += fv;
            });

            for (p, f) in particles.iter_mut().zip(forces.into_iter()) {
                p.force = f;
            }
            let electrostatic_energy =
                add_electrostatic_forces_particles(particles, box_length, &pme);

            //simulation_box.store_atoms_in_cells_particles(particles, &mut subcells, 10);

            // 4) recompute forces (LJ)
            //compute_forces_particles(particles, box_length, &mut subcells);

            // 5) velocity update (Verlet - second half step)
            for p in particles.iter_mut() {
                let a_new = p.force / p.mass;
                // where we do computing of the verlet
                p.update_velocity_verlet(a_new, dt);
            }

            // 6) measure temperature
            let dof = 3 * particles.len();
            let _system_temperature = compute_temperature_particles(&particles, dof);

            //println!("T = {system_temperature:.4}");

            // 6) thermostat (currently:  Berendsen and andersen  supported here)
            if thermostat == "berendsen" {
                apply_thermostat_berendsen_particles(particles, 300.0, 0.1, dt);
            } else if thermostat == "andersen" {
                apply_andersen_collisions(particles, 300.0, 1.0, dt);
            } else if thermostat == "nose_hoover" {
                apply_thermostat_nose_hoover_particles(
                    particles,
                    300.0,
                    10.0,
                    dt,
                    &mut xi_nose_hoover,
                )
            }

            // 7) recompute energy
            kinetic_energy = 0.0;
            for p in particles.iter() {
                kinetic_energy += 0.5 * p.mass * p.velocity.norm_squared();
            }
            potential_energy =
                site_site_energy_calculation(particles, box_length) + electrostatic_energy;
            total_energy = kinetic_energy + potential_energy;

            values.push(total_energy as f32);
        }

        info!(
            "Init particle energy | E_kin={kinetic_energy:.6} E_pot={potential_energy:.6} E_tot={total_energy:.6}"
        );

        // Optional: your running-average helper
        compute_average_val(&mut values, 2, number_of_steps as u64);
    }

    fn single_particle_energy(particles: &[Particle], idx: usize, box_length: f64) -> f64 {
        let mut energy = 0.0;
        let sigma_i = particles[idx].lj_parameters.sigma;
        let epsilon_i = particles[idx].lj_parameters.epsilon;

        for (j, other) in particles.iter().enumerate() {
            if j == idx {
                continue;
            }

            let sigma_j = other.lj_parameters.sigma;
            let epsilon_j = other.lj_parameters.epsilon;
            let sigma = (sigma_i + sigma_j) / 2.0;
            let epsilon = (epsilon_i * epsilon_j).sqrt();
            let r_vec = other.position - particles[idx].position;
            let r = minimum_image_convention(r_vec, box_length).norm();
            energy += lennard_jones_potential(r, sigma, epsilon);
        }

        energy
    }

    pub fn run_monte_carlo_particles(
        particles: &mut Vec<Particle>,
        number_of_steps: i32,
        box_length: f64,
        temperature: f64,
    ) {
        let mut values: Vec<f32> = Vec::new();
        let mut rng = rand::rng();
        let beta = 1.0 / temperature;
        let max_displacement = 0.05 * box_length;
        let mut accepted_moves: usize = 0;
        let mut attempted_moves: usize = 0;

        for _step in 0..number_of_steps {
            for idx in 0..particles.len() {
                let previous_position = particles[idx].position;
                let previous_energy = single_particle_energy(particles, idx, box_length);

                let displacement = Vector3::new(
                    rng.random_range(-max_displacement..max_displacement),
                    rng.random_range(-max_displacement..max_displacement),
                    rng.random_range(-max_displacement..max_displacement),
                );

                particles[idx].position += displacement;
                particles[idx].position.x = particles[idx].position.x.rem_euclid(box_length);
                particles[idx].position.y = particles[idx].position.y.rem_euclid(box_length);
                particles[idx].position.z = particles[idx].position.z.rem_euclid(box_length);

                let trial_energy = single_particle_energy(particles, idx, box_length);
                let delta_energy = trial_energy - previous_energy;
                let metropolis = (-beta * delta_energy).exp();

                let accepted = delta_energy <= 0.0 || rng.random::<f64>() < metropolis;
                attempted_moves += 1;

                if accepted {
                    accepted_moves += 1;
                } else {
                    particles[idx].position = previous_position;
                }
            }

            let potential_energy = site_site_energy_calculation(particles, box_length);
            values.push(potential_energy as f32);
        }

        let acceptance = accepted_moves as f64 / attempted_moves.max(1) as f64;
        info!(
            "Monte Carlo complete | attempted={} accepted={} ratio={acceptance:.4}",
            attempted_moves, accepted_moves
        );
        compute_average_val(&mut values, 2, number_of_steps as u64);
    }

    #[cfg(feature = "mpi")]
    fn rank_bounds(len: usize, rank: i32, size: i32) -> (usize, usize) {
        let rank = rank as usize;
        let size = size as usize;
        let base = len / size;
        let rem = len % size;

        let start = rank * base + rank.min(rem);
        let count = base + usize::from(rank < rem);
        (start, start + count)
    }

    #[cfg(feature = "mpi")]
    fn sync_particle_positions_and_velocities<C>(particles: &mut [Particle], world: &C)
    where
        C: mpi::topology::Communicator,
    {
        let n = particles.len();
        let mut packed = vec![0.0_f64; n * 6];

        if world.rank() == 0 {
            for (idx, particle) in particles.iter().enumerate() {
                let offset = idx * 6;
                packed[offset] = particle.position.x;
                packed[offset + 1] = particle.position.y;
                packed[offset + 2] = particle.position.z;
                packed[offset + 3] = particle.velocity.x;
                packed[offset + 4] = particle.velocity.y;
                packed[offset + 5] = particle.velocity.z;
            }
        }

        world.process_at_rank(0).broadcast_into(&mut packed[..]);

        for (idx, particle) in particles.iter_mut().enumerate() {
            let offset = idx * 6;
            particle.position.x = packed[offset];
            particle.position.y = packed[offset + 1];
            particle.position.z = packed[offset + 2];
            particle.velocity.x = packed[offset + 3];
            particle.velocity.y = packed[offset + 4];
            particle.velocity.z = packed[offset + 5];
        }
    }

    #[cfg(feature = "mpi")]
    pub fn compute_forces_particles_mpi<C>(
        particles: &mut Vec<Particle>,
        box_length: f64,
        world: &C,
    ) -> f64
    where
        C: mpi::topology::Communicator + mpi::traits::CommunicatorCollectives,
    {
        let n = particles.len();
        let (start, end) = rank_bounds(n, world.rank(), world.size());

        let mut local_forces = vec![0.0_f64; n * 3];
        let mut global_forces = vec![0.0_f64; n * 3];
        let mut local_potential = 0.0_f64;
        let mut global_potential = 0.0_f64;

        for i in start..end {
            for j in (i + 1)..n {
                let r_vec = particles[i].position - particles[j].position;
                let r_mic = minimum_image_convention(r_vec, box_length);
                let r = r_mic.norm();
                if r <= 1e-12 {
                    continue;
                }

                let si = particles[i].lj_parameters.sigma;
                let ei = particles[i].lj_parameters.epsilon;
                let sj = particles[j].lj_parameters.sigma;
                let ej = particles[j].lj_parameters.epsilon;
                let sigma = 0.5 * (si + sj);
                let epsilon = (ei * ej).sqrt();

                let f_mag = lennard_jones_force_scalar(r, sigma, epsilon);
                let f_vec = (r_mic / r) * f_mag;

                let ioff = 3 * i;
                local_forces[ioff] -= f_vec.x;
                local_forces[ioff + 1] -= f_vec.y;
                local_forces[ioff + 2] -= f_vec.z;

                let joff = 3 * j;
                local_forces[joff] += f_vec.x;
                local_forces[joff + 1] += f_vec.y;
                local_forces[joff + 2] += f_vec.z;

                local_potential += lennard_jones_potential(r, sigma, epsilon);
            }
        }

        world.all_reduce_into(
            &local_forces[..],
            &mut global_forces[..],
            SystemOperation::sum(),
        );
        world.all_reduce_into(
            &local_potential,
            &mut global_potential,
            SystemOperation::sum(),
        );

        for (idx, particle) in particles.iter_mut().enumerate() {
            let offset = 3 * idx;
            particle.force = Vector3::new(
                global_forces[offset],
                global_forces[offset + 1],
                global_forces[offset + 2],
            );
        }

        global_potential
    }

    #[cfg(feature = "mpi")]
    pub fn run_md_nve_particles_mpi<C>(
        particles: &mut Vec<Particle>,
        number_of_steps: i32,
        dt: f64,
        box_length: f64,
        thermostat: &str,
        world: &C,
    ) where
        C: mpi::topology::Communicator + mpi::traits::CommunicatorCollectives,
    {
        sync_particle_positions_and_velocities(particles, world);

        let mut values: Vec<f32> = Vec::new();
        let mut potential_energy = compute_forces_particles_mpi(particles, box_length, world);

        let n = particles.len();
        let (start, end) = rank_bounds(n, world.rank(), world.size());
        let local_kinetic: f64 = particles[start..end]
            .iter()
            .map(|p| 0.5 * p.mass * p.velocity.norm_squared())
            .sum();
        let mut kinetic_energy = 0.0;
        world.all_reduce_into(&local_kinetic, &mut kinetic_energy, SystemOperation::sum());
        let mut total_energy = kinetic_energy + potential_energy;

        if world.rank() == 0 {
            info!(
                "Init particle energy (MPI) | E_kin={kinetic_energy:.6} E_pot={potential_energy:.6} E_tot={total_energy:.6}"
            );
        }

        let mut xi_nose_hoover = 0.0;

        for _step in 0..number_of_steps {
            let mut a_old: Vec<Vector3<f64>> = Vec::with_capacity(particles.len());
            for p in particles.iter() {
                a_old.push(p.force / p.mass);
            }

            for (atom, a_o) in particles.iter_mut().zip(a_old.iter()) {
                atom.velocity += 0.5 * a_o * dt;
            }

            for atom in particles.iter_mut() {
                atom.update_position_verlet(dt);
            }

            pbc_update(particles, box_length);
            potential_energy = compute_forces_particles_mpi(particles, box_length, world);

            for p in particles.iter_mut() {
                let a_new = p.force / p.mass;
                p.update_velocity_verlet(a_new, dt);
            }

            if thermostat == "berendsen" {
                apply_thermostat_berendsen_particles(particles, 300.0, 0.1, dt);
            } else if thermostat == "andersen" {
                apply_andersen_collisions(particles, 300.0, 1.0, dt);
            } else if thermostat == "nose_hoover" {
                apply_thermostat_nose_hoover_particles(
                    particles,
                    300.0,
                    10.0,
                    dt,
                    &mut xi_nose_hoover,
                );
            }

            // Keep all ranks synchronized even when stochastic thermostats are used.
            sync_particle_positions_and_velocities(particles, world);

            let local_kinetic: f64 = particles[start..end]
                .iter()
                .map(|p| 0.5 * p.mass * p.velocity.norm_squared())
                .sum();
            world.all_reduce_into(&local_kinetic, &mut kinetic_energy, SystemOperation::sum());
            total_energy = kinetic_energy + potential_energy;

            if world.rank() == 0 {
                values.push(total_energy as f32);
            }
        }

        if world.rank() == 0 {
            compute_average_val(&mut values, 2, number_of_steps as u64);
        }
    }

    /*

    Systems portion of the code

     */

    pub fn run_md_nve_systems(
        systems: &mut Vec<System>,
        number_of_steps: i32,
        dt: f64,
        box_length: f64,
        thermostat: &str,
    ) {
        let mut values: Vec<f32> = Vec::new();
        let mut total_energy = 0.0;
        let mut kinetic_energy = 0.0;
        let mut potential_energy = 0.0;
        let pme = PmeConfig::default();

        // Create the subcells for the simulation box
        let simulation_box = cell_subdivision::SimulationBox {
            x_dimension: box_length,
            y_dimension: box_length,
            z_dimension: box_length,
        };
        // Create the subcells - here we have used a subdivision of 10 for the cells
        let mut subcells = simulation_box.create_subcells(10);
        // Store the coordinates in cells
        simulation_box.store_atoms_in_cells_systems(systems, &mut subcells, 10);
        // --- initial forces and energy ---

        info!(
            "Init systems energy | E_kin={kinetic_energy:.6} E_pot={potential_energy:.6} E_tot={total_energy:.6}"
        );

        for sys in systems.iter_mut() {
            for a in sys.atoms.iter_mut() {
                a.force = Vector3::zeros();
            }
            compute_bonded_forces(
                &mut sys.atoms,
                &sys.bonds,
                &sys.angles,
                &sys.dihedrals,
                &sys.impropers,
                box_length,
            );
        }
        compute_intermolecular_forces_systems(systems, box_length);
        let _ = add_electrostatic_forces_systems(systems, box_length, &pme);

        // this is only used if we apply nose hoover
        let mut xi_nose_hoover = vec![0.0; systems.len()];

        // --- time integration loop ---
        for _step in 0..number_of_steps {
            for sys in systems.iter_mut() {
                let mut a_old: Vec<Vector3<f64>> = Vec::with_capacity(sys.atoms.len());

                for a in sys.atoms.iter() {
                    a_old.push(a.force / a.mass);
                }

                for (atom, a_o) in sys.atoms.iter_mut().zip(a_old.iter()) {
                    atom.velocity += 0.5 * a_o * dt;
                }

                for atom in sys.atoms.iter_mut() {
                    atom.update_position_verlet(dt);
                }

                pbc_update(&mut sys.atoms, box_length);

                for a in sys.atoms.iter_mut() {
                    a.force = Vector3::zeros();
                }
                compute_bonded_forces(
                    &mut sys.atoms,
                    &sys.bonds,
                    &sys.angles,
                    &sys.dihedrals,
                    &sys.impropers,
                    box_length,
                );
            }

            compute_intermolecular_forces_systems(systems, box_length);
            let electrostatic_energy = add_electrostatic_forces_systems(systems, box_length, &pme);

            for (s, sys) in systems.iter_mut().enumerate() {
                for a in sys.atoms.iter_mut() {
                    let a_new = a.force / a.mass;
                    a.update_velocity_verlet(a_new, dt);
                }

                let dof = 3 * sys.atoms.len();
                let _system_temperature = compute_temperature_particles(&sys.atoms, dof);
                if thermostat == "berendsen" {
                    apply_thermostat_berendsen_particles(&mut sys.atoms, 300.0, 0.1, dt);
                } else if thermostat == "nose_hoover" {
                    apply_thermostat_nose_hoover_particles(
                        &mut sys.atoms,
                        300.0,
                        10.0,
                        dt,
                        &mut xi_nose_hoover[s],
                    )
                }
            }

            kinetic_energy = 0.0;
            potential_energy = 0.0;
            for sys in systems.iter_mut() {
                for a in sys.atoms.iter() {
                    kinetic_energy += 0.5 * a.mass * a.velocity.norm_squared();
                }
                potential_energy += compute_bonded_energy(
                    &mut sys.atoms,
                    &sys.bonds,
                    &sys.angles,
                    &sys.dihedrals,
                    &sys.impropers,
                    box_length,
                );
            }
            potential_energy += intermolecular_site_site_energy_systems(systems, box_length);
            potential_energy += electrostatic_energy;

            total_energy = kinetic_energy + potential_energy;
            values.push(total_energy as f32);
            info!("Step {_step:>4} | E_tot={total_energy:.6} E_kin={kinetic_energy:.6} E_pot={potential_energy:.6}");
        }

        compute_average_val(&mut values, 2, number_of_steps as u64);
    }

    pub fn run_md_nve(
        state: &mut InitOutput,
        number_of_steps: i32,
        dt: f64,
        box_length: f64,
        thermostat: &str,
        cutoff: f64,
    ) {
        if thermostat == "monte_carlo" {
            match state {
                InitOutput::Particles(particles) => {
                    run_monte_carlo_particles(particles, number_of_steps, box_length, 300.0);
                }
                InitOutput::Systems(_systems) => {
                    info!(
                        "Monte Carlo mode is currently supported for particle simulations only; using velocity Verlet for systems."
                    );
                }
            }
        }

        match state {
            InitOutput::Particles(particles) => {
                if thermostat == "monte_carlo" {
                    return;
                }
                run_md_nve_particles(particles, number_of_steps, dt, box_length, thermostat, 30.0);
            }
            InitOutput::Systems(systems) => {
                run_md_nve_systems(systems, number_of_steps, dt, box_length, thermostat);
            }
        }
    }

    #[cfg(feature = "mpi")]
    pub fn run_md_nve_mpi<C>(
        state: &mut InitOutput,
        number_of_steps: i32,
        dt: f64,
        box_length: f64,
        thermostat: &str,
        world: &C,
    ) where
        C: mpi::topology::Communicator + mpi::traits::CommunicatorCollectives,
    {
        match state {
            InitOutput::Particles(particles) => {
                run_md_nve_particles_mpi(
                    particles,
                    number_of_steps,
                    dt,
                    box_length,
                    thermostat,
                    world,
                );
            }
            InitOutput::Systems(systems) => {
                if world.rank() == 0 {
                    info!(
                        "MPI NVE currently supports particle systems; falling back to serial systems integration."
                    );
                }
                run_md_nve_systems(systems, number_of_steps, dt, box_length, thermostat);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use log::{error, info};

    // lennard-jones double loop test
    #[test]
    fn test_double_loop() {
        let _lj_params = lennard_jones_simulations::LJParameters {
            epsilon: 1.0,
            sigma: 1.0,
            number_of_atoms: 2,
        };
        // call the double_loop function

        // assert that the result is as expected
        //assert_eq!(result, expected_result)
    }

    #[test]
    fn test_lennard_jones() {
        let _sigma = 1.0;
        let _epsilon = 1.0;
        let _lj_params_new = lennard_jones_simulations::LJParameters {
            epsilon: 1.0,
            sigma: 1.0,
            number_of_atoms: 10,
        };
    }

    #[test]
    fn berenden_pull_towards_target() {
        /* mock velocities - T = 300K
           let mut t = 300.0;
           let t0 = 350.0;
           let dt = 0.001;
           let tau = 0.1;
        */
        let t0 = 300.0;

        // Define the new simulation for nve
        let mut new_simulation_md =
            match lennard_jones_simulations::create_atoms_with_set_positions_and_velocities(
                10, 300.0, 30.0, 10.0, 10.0, false,
            ) {
                // How to handle errors - we are returning a result or a string
                Ok(atoms) => atoms,
                Err(e) => {
                    error!("Failed to create atoms: {e}");
                    return; // Exit early or handle the error as needed
                }
            };
        lennard_jones_simulations::run_md_nve(
            &mut new_simulation_md,
            1000,
            0.5,
            10.0,
            "berendsen",
            30.0,
        );

        let dof = match &new_simulation_md {
            lennard_jones_simulations::InitOutput::Particles(particles) => 3 * particles.len(),
            lennard_jones_simulations::InitOutput::Systems(systems) => {
                3 * systems.iter().map(|sys| sys.atoms.len()).sum::<usize>()
            }
        };
        //let dof = 3 * new_simulation_md.len();
        // compute the final temperature of the system
        let t = lennard_jones_simulations::compute_temperature(&mut new_simulation_md, dof);
        info!("Final temperature={t:.3}, target={t0:.3}");
        assert!((t - t0).abs() < 5.0, "Temperature should approach target");
    }
}
