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

/// Coulomb prefactor for MD units used across FerrumMD.
///
/// Unit convention:
/// - distance: nm
/// - charge: elementary charge (e)
/// - energy: kJ/mol
///
/// The prefactor is k_e = 138.935457644 kJ mol^-1 nm e^-2, which is the
/// standard molecular-simulation conversion for 1/(4*pi*epsilon_0).
#[inline]
fn coulomb_prefactor() -> f64 {
    138.935_457_644
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
            real_cutoff: 1.2,
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
    use nalgebra::{zero, Matrix3, Vector3};
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
    use crate::molecule::shake_rattle::shake_rattle::{self, DistanceConstraint};

    use crate::lennard_jones_simulations::cell_subdivision::MolecularCoordinates;

    #[derive(Clone, Debug)]
    pub struct LJParameters {
        // lennard jones parameters and the number of atoms that we have of that parameter
        pub epsilon: f64,
        pub sigma: f64,
        pub number_of_atoms: i32,
    }

    #[derive(Clone, Debug)]
    /// Particle state in explicit MD units.
    ///
    /// - `position`: nm
    /// - `velocity`: nm/ps
    /// - `force`: kJ/(mol·nm)
    /// - `mass`: amu
    /// - `charge`: elementary charge (e)
    /// - `energy`: kJ/mol
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

    #[derive(Clone, Copy, Debug)]
    /// Configuration for system-level simulation loops.
    ///
    /// - `cutoff`, `neighbor_skin`: nm
    /// - `pme.alpha`: nm^-1
    /// - integration timestep (`dt` arguments in run functions): ps
    pub struct SystemSimulationConfig {
        pub cutoff: f64,
        pub neighbor_skin: f64,
        pub neighbor_rebuild_interval: usize,
        pub pme: PmeConfig,
    }

    impl Default for SystemSimulationConfig {
        fn default() -> Self {
            Self {
                cutoff: 1.2,
                neighbor_skin: 0.2,
                neighbor_rebuild_interval: 10,
                pme: PmeConfig::default(),
            }
        }
    }

    #[derive(Clone, Debug)]
    struct SystemVerletNeighborList {
        pairs: Vec<(usize, usize)>,
        reference_positions: Vec<Vector3<f64>>,
        cutoff: f64,
        skin: f64,
        rebuild_interval: usize,
        steps_since_rebuild: usize,
    }

    impl SystemVerletNeighborList {
        fn new(cutoff: f64, skin: f64, rebuild_interval: usize) -> Self {
            Self {
                pairs: Vec::new(),
                reference_positions: Vec::new(),
                cutoff,
                skin,
                rebuild_interval: rebuild_interval.max(1),
                steps_since_rebuild: 0,
            }
        }

        fn capture_positions(systems: &[System]) -> Vec<Vector3<f64>> {
            systems
                .iter()
                .flat_map(|sys| sys.atoms.iter().map(|a| a.position))
                .collect()
        }

        fn atom_count(systems: &[System]) -> usize {
            systems.iter().map(|sys| sys.atoms.len()).sum()
        }

        fn max_displacement_since_rebuild(&self, systems: &[System], box_length: f64) -> f64 {
            if self.reference_positions.is_empty() {
                return f64::INFINITY;
            }

            let mut max_disp: f64 = 0.0;
            let mut idx = 0usize;
            for sys in systems.iter() {
                for atom in sys.atoms.iter() {
                    if idx >= self.reference_positions.len() {
                        return f64::INFINITY;
                    }
                    let dr = minimum_image_convention(
                        atom.position - self.reference_positions[idx],
                        box_length,
                    );
                    max_disp = max_disp.max(dr.norm());
                    idx += 1;
                }
            }

            max_disp
        }

        fn rebuild_if_needed(&mut self, systems: &[System], box_length: f64) {
            let atom_count = Self::atom_count(systems);
            let needs_size_rebuild = self.reference_positions.len() != atom_count;
            let max_disp = self.max_displacement_since_rebuild(systems, box_length);
            let needs_skin_rebuild = max_disp >= 0.5 * self.skin;
            let needs_interval_rebuild = self.steps_since_rebuild >= self.rebuild_interval;

            if needs_size_rebuild || needs_skin_rebuild || needs_interval_rebuild {
                self.rebuild_pairs(systems, box_length);
            } else {
                self.steps_since_rebuild += 1;
            }
        }

        fn rebuild_pairs(&mut self, systems: &[System], box_length: f64) {
            self.pairs.clear();
            self.reference_positions = Self::capture_positions(systems);
            self.steps_since_rebuild = 0;

            let mut atom_map: Vec<(usize, usize)> =
                Vec::with_capacity(self.reference_positions.len());
            for (sys_idx, sys) in systems.iter().enumerate() {
                for atom_idx in 0..sys.atoms.len() {
                    atom_map.push((sys_idx, atom_idx));
                }
            }

            let neighbor_cutoff = self.cutoff + self.skin;
            let neighbor_cutoff2 = neighbor_cutoff * neighbor_cutoff;

            for i in 0..atom_map.len() {
                for j in (i + 1)..atom_map.len() {
                    if atom_map[i].0 == atom_map[j].0 {
                        continue; // omit intramolecular pairs
                    }

                    let ri = self.reference_positions[i];
                    let rj = self.reference_positions[j];
                    let rij = minimum_image_convention(rj - ri, box_length);
                    if rij.norm_squared() <= neighbor_cutoff2 {
                        self.pairs.push((i, j));
                    }
                }
            }
        }
    }

    pub enum InitOutput {
        Particles(Vec<Particle>), // define a particles system (single point particle)
        Systems(Vec<System>),     // Define actual molecules
    }

    pub enum InitMode {
        Atoms,
        Molecules,
    }

    #[derive(Clone, Debug)]
    pub enum ConstraintMode {
        Flexible,
        ShakeRattle,
        SettlePreferred,
    }

    #[derive(Clone, Debug)]
    pub struct ConstraintOptions {
        pub mode: ConstraintMode,
        pub constraints_by_system: Vec<Vec<DistanceConstraint>>,
        pub tolerance: f64,
        pub max_iter: usize,
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
        const MIN_SEPARATION: f64 = 0.35;
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
                let mut used_translations: Vec<Vector3<f64>> = Vec::new();
                for sys in systems.iter_mut() {
                    // Randomize each molecule position as a rigid translation so cloned
                    // systems do not remain perfectly overlapped at initialization.
                    let translation = sample_position_with_min_separation(
                        &mut rng,
                        box_length,
                        MIN_SEPARATION,
                        &used_translations,
                    );
                    used_translations.push(translation);

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
        let mut existing_positions: Vec<Vector3<f64>> = Vec::new();
        let mut rng = rand::rng();
        // Create the number of atoms in the system with the system as necessary
        if !use_atom {
            for index in 0..number_of_atoms {
                // For each particle, we wish to create the initial posiiton,
                // the velocity, and the LJ parameters attached to it
                let mut particle = Particle {
                    // create position for the atom in question
                    position: sample_position_with_min_separation(
                        &mut rng,
                        box_dim_max,
                        0.9,
                        &existing_positions,
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
                existing_positions.push(particle.position);
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

    fn sample_position_with_min_separation(
        rng: &mut impl Rng,
        box_length: f64,
        min_separation: f64,
        existing_positions: &[Vector3<f64>],
    ) -> Vector3<f64> {
        let min_sep2 = min_separation * min_separation;
        const MAX_TRIES: usize = 10_000;
        for _ in 0..MAX_TRIES {
            let candidate = Vector3::new(
                rng.random_range(0.0..box_length),
                rng.random_range(0.0..box_length),
                rng.random_range(0.0..box_length),
            );
            let has_overlap = existing_positions.iter().any(|other| {
                let dr = minimum_image_convention(candidate - *other, box_length);
                dr.norm_squared() < min_sep2
            });
            if !has_overlap {
                return candidate;
            }
        }
        Vector3::new(
            rng.random_range(0.0..box_length),
            rng.random_range(0.0..box_length),
            rng.random_range(0.0..box_length),
        )
    }

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

    fn flatten_system_atom_map(systems: &[System]) -> Vec<(usize, usize)> {
        let mut atom_map = Vec::new();
        for (sys_idx, sys) in systems.iter().enumerate() {
            for atom_idx in 0..sys.atoms.len() {
                atom_map.push((sys_idx, atom_idx));
            }
        }
        atom_map
    }

    fn lj_coulomb_scaling_for_pair(
        systems: &[System],
        atom_map: &[(usize, usize)],
        global_i: usize,
        global_j: usize,
    ) -> Option<(f64, f64)> {
        let (sys_i, atom_i) = atom_map[global_i];
        let (sys_j, atom_j) = atom_map[global_j];
        if sys_i != sys_j {
            return Some((1.0, 1.0));
        }
        systems[sys_i]
            .nonbonded_scaling_for_pair(atom_i, atom_j)
            .map(|scale| (scale.lj_scale, scale.coulomb_scale))
    }

    fn compute_intermolecular_forces_systems_with_pairs(
        systems: &mut [System],
        box_length: f64,
        pairs: &[(usize, usize)],
        cutoff: f64,
    ) -> f64 {
        let mut total_energy = 0.0;
        let atom_map = flatten_system_atom_map(systems);
        let cutoff2 = cutoff * cutoff;

        for &(global_i, global_j) in pairs.iter() {
            let (sys_i_idx_raw, atom_i_idx_raw) = atom_map[global_i];
            let (sys_j_idx_raw, atom_j_idx_raw) = atom_map[global_j];
            if sys_i_idx_raw == sys_j_idx_raw {
                continue;
            }
            let Some((lj_scale, _)) =
                lj_coulomb_scaling_for_pair(systems, &atom_map, global_i, global_j)
            else {
                continue;
            };

            let (sys_i_idx, atom_i_idx, sys_j_idx, atom_j_idx, swap_sign) =
                if sys_i_idx_raw < sys_j_idx_raw {
                    (
                        sys_i_idx_raw,
                        atom_i_idx_raw,
                        sys_j_idx_raw,
                        atom_j_idx_raw,
                        1.0,
                    )
                } else {
                    (
                        sys_j_idx_raw,
                        atom_j_idx_raw,
                        sys_i_idx_raw,
                        atom_i_idx_raw,
                        -1.0,
                    )
                };

            let (left, right) = systems.split_at_mut(sys_j_idx);
            let sys_i = &mut left[sys_i_idx];
            let sys_j = &mut right[0];
            let atom_i = &mut sys_i.atoms[atom_i_idx];
            let atom_j = &mut sys_j.atoms[atom_j_idx];

            let r_vec = atom_j.position - atom_i.position;
            let r_mic = minimum_image_convention(r_vec, box_length);
            let r2 = r_mic.norm_squared();
            if r2 <= 1e-24 || r2 > cutoff2 {
                continue;
            }
            let r = safe_norm(r2.sqrt());

            let sigma = 0.5 * (atom_i.lj_parameters.sigma + atom_j.lj_parameters.sigma);
            let epsilon = (atom_i.lj_parameters.epsilon * atom_j.lj_parameters.epsilon).sqrt();

            let f_mag = lennard_jones_force_scalar(r, sigma, epsilon);
            let f_vec = (r_mic / r) * f_mag * swap_sign;

            atom_i.force -= f_vec * lj_scale;
            atom_j.force += f_vec * lj_scale;

            total_energy += lj_scale * lennard_jones_potential(r, sigma, epsilon);
        }

        total_energy +=
            compute_intramolecular_pair_lj_forces_and_energy(systems, box_length, cutoff);

        total_energy
    }

    fn compute_intramolecular_pair_lj_forces_and_energy(
        systems: &mut [System],
        box_length: f64,
        cutoff: f64,
    ) -> f64 {
        let cutoff2 = cutoff * cutoff;
        let mut energy = 0.0;

        for sys in systems.iter_mut() {
            let pairs: Vec<(
                (usize, usize),
                crate::molecule::molecule::NonbondedPairScaling,
            )> = sys.pair_scalings.iter().map(|(k, v)| (*k, *v)).collect();
            for ((i, j), scaling) in pairs {
                if scaling.lj_scale == 0.0 || sys.exclusions.contains(&(i, j)) {
                    continue;
                }
                if i >= sys.atoms.len() || j >= sys.atoms.len() || i == j {
                    continue;
                }
                let (ai, aj) = if i < j {
                    let (left, right) = sys.atoms.split_at_mut(j);
                    (&mut left[i], &mut right[0])
                } else {
                    let (left, right) = sys.atoms.split_at_mut(i);
                    (&mut right[0], &mut left[j])
                };
                let r_vec = aj.position - ai.position;
                let r_mic = minimum_image_convention(r_vec, box_length);
                let r2 = r_mic.norm_squared();
                if r2 <= 1e-24 || r2 > cutoff2 {
                    continue;
                }
                let r = safe_norm(r2.sqrt());
                let sigma = 0.5 * (ai.lj_parameters.sigma + aj.lj_parameters.sigma);
                let epsilon = (ai.lj_parameters.epsilon * aj.lj_parameters.epsilon).sqrt();
                let f_mag = lennard_jones_force_scalar(r, sigma, epsilon);
                let f_vec = (r_mic / r) * f_mag * scaling.lj_scale;
                ai.force -= f_vec;
                aj.force += f_vec;
                energy += scaling.lj_scale * lennard_jones_potential(r, sigma, epsilon);
            }
        }

        energy
    }

    pub fn compute_intermolecular_forces_systems(systems: &mut [System], box_length: f64) -> f64 {
        let cutoff = 0.5 * box_length;
        compute_intermolecular_forces_systems_cutoff(systems, box_length, cutoff)
    }

    pub fn compute_intermolecular_forces_systems_cutoff(
        systems: &mut [System],
        box_length: f64,
        cutoff: f64,
    ) -> f64 {
        /*
        Compute Lennard-Jones interactions between atoms belonging to different systems.
        Intra-molecular interactions are omitted here and handled by bonded terms.
         */
        let mut total_energy = 0.0;
        let cutoff2 = cutoff * cutoff;

        for i in 0..systems.len() {
            for j in (i + 1)..systems.len() {
                let (left, right) = systems.split_at_mut(j);
                let sys_i = &mut left[i];
                let sys_j = &mut right[0];

                for atom_i in sys_i.atoms.iter_mut() {
                    for atom_j in sys_j.atoms.iter_mut() {
                        let r_vec = atom_j.position - atom_i.position;
                        let r_mic = minimum_image_convention(r_vec, box_length);
                        let r2 = r_mic.norm_squared();
                        if r2 > cutoff2 {
                            continue;
                        }
                        let r = safe_norm(r2.sqrt());

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
        let cutoff = 0.5 * box_length;
        intermolecular_site_site_energy_systems_cutoff(systems, box_length, cutoff)
    }

    pub fn intermolecular_site_site_energy_systems_cutoff(
        systems: &[System],
        box_length: f64,
        cutoff: f64,
    ) -> f64 {
        /*
        Compute Lennard-Jones potential energy between atoms in different systems.
         */
        let mut total_energy = 0.0;
        let cutoff2 = cutoff * cutoff;

        for i in 0..systems.len() {
            for j in (i + 1)..systems.len() {
                let sys_i = &systems[i];
                let sys_j = &systems[j];

                for atom_i in sys_i.atoms.iter() {
                    for atom_j in sys_j.atoms.iter() {
                        let r_vec = atom_j.position - atom_i.position;
                        let r_mic = minimum_image_convention(r_vec, box_length);
                        let r2 = r_mic.norm_squared();
                        if r2 > cutoff2 {
                            continue;
                        }
                        let r = safe_norm(r2.sqrt());

                        let sigma = 0.5 * (atom_i.lj_parameters.sigma + atom_j.lj_parameters.sigma);
                        let epsilon =
                            (atom_i.lj_parameters.epsilon * atom_j.lj_parameters.epsilon).sqrt();

                        total_energy += lennard_jones_potential(r, sigma, epsilon);
                    }
                }
            }
        }

        total_energy += intramolecular_pair_lj_energy_systems(systems, box_length, cutoff);

        total_energy
    }

    fn intramolecular_pair_lj_energy_systems(
        systems: &[System],
        box_length: f64,
        cutoff: f64,
    ) -> f64 {
        let cutoff2 = cutoff * cutoff;
        let mut energy = 0.0;
        for sys in systems {
            for (&(i, j), scaling) in &sys.pair_scalings {
                if scaling.lj_scale == 0.0 || sys.exclusions.contains(&(i, j)) {
                    continue;
                }
                if i >= sys.atoms.len() || j >= sys.atoms.len() || i == j {
                    continue;
                }
                let r_vec = sys.atoms[j].position - sys.atoms[i].position;
                let r_mic = minimum_image_convention(r_vec, box_length);
                let r2 = r_mic.norm_squared();
                if r2 <= 1e-24 || r2 > cutoff2 {
                    continue;
                }
                let r = safe_norm(r2.sqrt());
                let sigma =
                    0.5 * (sys.atoms[i].lj_parameters.sigma + sys.atoms[j].lj_parameters.sigma);
                let epsilon = (sys.atoms[i].lj_parameters.epsilon
                    * sys.atoms[j].lj_parameters.epsilon)
                    .sqrt();
                energy += scaling.lj_scale * lennard_jones_potential(r, sigma, epsilon);
            }
        }
        energy
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

    pub fn compute_electrostatic_forces_systems_with_config(
        systems: &mut [System],
        box_length: f64,
        pme: &PmeConfig,
    ) -> f64 {
        add_electrostatic_forces_systems(systems, box_length, pme)
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
        let atom_map = flatten_system_atom_map(systems);
        let mut energy = 0.0;
        let alpha = pme.alpha;
        let rc = pme.real_cutoff;
        let k_e = coulomb_prefactor();

        for global_i in 0..atom_map.len() {
            for global_j in (global_i + 1)..atom_map.len() {
                let Some((_, coulomb_scale)) =
                    lj_coulomb_scaling_for_pair(systems, &atom_map, global_i, global_j)
                else {
                    continue;
                };
                if coulomb_scale == 0.0 {
                    continue;
                }

                let (sys_i_idx_raw, atom_i_idx_raw) = atom_map[global_i];
                let (sys_j_idx_raw, atom_j_idx_raw) = atom_map[global_j];

                let (sys_i_idx, atom_i_idx, sys_j_idx, atom_j_idx, swap_sign) =
                    if sys_i_idx_raw < sys_j_idx_raw {
                        (
                            sys_i_idx_raw,
                            atom_i_idx_raw,
                            sys_j_idx_raw,
                            atom_j_idx_raw,
                            1.0,
                        )
                    } else if sys_j_idx_raw < sys_i_idx_raw {
                        (
                            sys_j_idx_raw,
                            atom_j_idx_raw,
                            sys_i_idx_raw,
                            atom_i_idx_raw,
                            -1.0,
                        )
                    } else if atom_i_idx_raw < atom_j_idx_raw {
                        (
                            sys_i_idx_raw,
                            atom_i_idx_raw,
                            sys_j_idx_raw,
                            atom_j_idx_raw,
                            1.0,
                        )
                    } else {
                        (
                            sys_j_idx_raw,
                            atom_j_idx_raw,
                            sys_i_idx_raw,
                            atom_i_idx_raw,
                            -1.0,
                        )
                    };

                if sys_i_idx == sys_j_idx {
                    let sys = &mut systems[sys_i_idx];
                    let (atom_i, atom_j) = if atom_i_idx < atom_j_idx {
                        let (left, right) = sys.atoms.split_at_mut(atom_j_idx);
                        (&mut left[atom_i_idx], &mut right[0])
                    } else {
                        let (left, right) = sys.atoms.split_at_mut(atom_i_idx);
                        (&mut right[0], &mut left[atom_j_idx])
                    };
                    let qi = atom_i.charge;
                    let qj = atom_j.charge;
                    if qi == 0.0 && qj == 0.0 {
                        continue;
                    }
                    let rij =
                        minimum_image_convention(atom_j.position - atom_i.position, box_length);
                    let r = rij.norm();
                    if r <= 1e-12 || r > rc {
                        continue;
                    }
                    let ar = alpha * r;
                    let erfc_ar = erfc_approx(ar);
                    let exp_term = (-(ar * ar)).exp();
                    let qq = k_e * qi * qj * coulomb_scale;
                    energy += qq * erfc_ar / r;
                    let scalar = qq
                        * (erfc_ar / (r * r)
                            + (2.0 * alpha / std::f64::consts::PI.sqrt()) * exp_term / r);
                    let f_vec = -(rij / r) * scalar * swap_sign;
                    atom_i.force += f_vec;
                    atom_j.force -= f_vec;
                    continue;
                }

                let (left, right) = systems.split_at_mut(sys_j_idx);
                let sys_i = &mut left[sys_i_idx];
                let sys_j = &mut right[0];
                let atom_i = &mut sys_i.atoms[atom_i_idx];
                let atom_j = &mut sys_j.atoms[atom_j_idx];
                let qi = atom_i.charge;
                let qj = atom_j.charge;
                if qi == 0.0 && qj == 0.0 {
                    continue;
                }

                let rij = minimum_image_convention(atom_j.position - atom_i.position, box_length);
                let r = rij.norm();
                if r <= 1e-12 || r > rc {
                    continue;
                }

                let ar = alpha * r;
                let erfc_ar = erfc_approx(ar);
                let exp_term = (-(ar * ar)).exp();
                let qq = k_e * qi * qj * coulomb_scale;
                energy += qq * erfc_ar / r;

                let scalar = qq
                    * (erfc_ar / (r * r)
                        + (2.0 * alpha / std::f64::consts::PI.sqrt()) * exp_term / r);
                let f_vec = -(rij / r) * scalar * swap_sign;
                atom_i.force += f_vec;
                atom_j.force -= f_vec;
            }
        }

        energy + add_electrostatic_reciprocal_systems(systems, box_length, pme, &atom_map)
    }

    fn add_electrostatic_reciprocal_systems(
        systems: &mut [System],
        box_length: f64,
        pme: &PmeConfig,
        atom_map: &[(usize, usize)],
    ) -> f64 {
        if atom_map.is_empty() {
            return 0.0;
        }

        let volume = box_length.powi(3);
        if volume <= 0.0 {
            return 0.0;
        }

        let alpha = pme.alpha;
        let kmax = pme.kmax;
        let k_e = coulomb_prefactor();
        let two_pi_over_l = 2.0 * std::f64::consts::PI / box_length;

        let positions: Vec<Vector3<f64>> = atom_map
            .iter()
            .map(|&(sys_idx, atom_idx)| systems[sys_idx].atoms[atom_idx].position)
            .collect();
        let charges: Vec<f64> = atom_map
            .iter()
            .map(|&(sys_idx, atom_idx)| systems[sys_idx].atoms[atom_idx].charge)
            .collect();

        let mut force_updates = vec![Vector3::<f64>::zeros(); atom_map.len()];
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
                    let coef_energy = (2.0 * std::f64::consts::PI * k_e / volume) * damp / k2;
                    let coef_force = (4.0 * std::f64::consts::PI * k_e / volume) * damp / k2;

                    for i in 0..atom_map.len() {
                        if charges[i] == 0.0 {
                            continue;
                        }
                        for j in (i + 1)..atom_map.len() {
                            if charges[j] == 0.0 {
                                continue;
                            }
                            let Some((_, coulomb_scale)) =
                                lj_coulomb_scaling_for_pair(systems, atom_map, i, j)
                            else {
                                continue;
                            };
                            if coulomb_scale == 0.0 {
                                continue;
                            }

                            let rij = positions[j] - positions[i];
                            let phase = kvec.dot(&rij);
                            let pair_charge = charges[i] * charges[j] * coulomb_scale;
                            energy += 2.0 * coef_energy * pair_charge * phase.cos();

                            let force_vec = coef_force * pair_charge * phase.sin() * kvec;
                            force_updates[i] -= force_vec;
                            force_updates[j] += force_vec;
                        }
                    }
                }
            }
        }

        let self_energy: f64 = charges.iter().map(|q| q * q).sum::<f64>()
            * (-k_e * alpha / std::f64::consts::PI.sqrt());
        energy += self_energy;

        for (global_idx, &(sys_idx, atom_idx)) in atom_map.iter().enumerate() {
            systems[sys_idx].atoms[atom_idx].force += force_updates[global_idx];
        }

        energy
    }

    fn enforce_settle_position(system: &mut System) {
        let Some(rigid) = system.rigid_water else {
            return;
        };
        let o = rigid.oxygen_index;
        let h1 = rigid.hydrogen1_index;
        let h2 = rigid.hydrogen2_index;
        if o >= system.atoms.len() || h1 >= system.atoms.len() || h2 >= system.atoms.len() {
            return;
        }

        let m_o = system.atoms[o].mass;
        let m_h1 = system.atoms[h1].mass;
        let m_h2 = system.atoms[h2].mass;
        let total_mass = m_o + m_h1 + m_h2;
        if total_mass <= 0.0 {
            return;
        }
        let com = (system.atoms[o].position * m_o
            + system.atoms[h1].position * m_h1
            + system.atoms[h2].position * m_h2)
            / total_mass;

        let mut e1 = system.atoms[h1].position - system.atoms[o].position;
        if e1.norm_squared() < 1e-18 {
            e1 = Vector3::new(1.0, 0.0, 0.0);
        }
        e1 = e1.normalize();
        let mut in_plane = system.atoms[h2].position - system.atoms[o].position;
        in_plane -= e1 * in_plane.dot(&e1);
        if in_plane.norm_squared() < 1e-18 {
            in_plane = Vector3::new(0.0, 1.0, 0.0);
        }
        let e2 = in_plane.normalize();
        let e3 = e1.cross(&e2);

        let d_oh = rigid.oh_distance;
        let d_hh = rigid.hh_distance;
        let x = (d_hh * d_hh) / (2.0 * d_oh);
        let y_sq = (d_oh * d_oh - x * x).max(0.0);
        let y = y_sq.sqrt();
        let r_h1 = Vector3::new(x, y, 0.0);
        let r_h2 = Vector3::new(x, -y, 0.0);
        let r_o = -(m_h1 * r_h1 + m_h2 * r_h2) / m_o;

        let world = |r: Vector3<f64>| com + e1 * r.x + e2 * r.y + e3 * r.z;
        system.atoms[o].position = world(r_o);
        system.atoms[h1].position = world(r_h1);
        system.atoms[h2].position = world(r_h2);
    }

    fn enforce_settle_velocity(system: &mut System) {
        let Some(rigid) = system.rigid_water else {
            return;
        };
        let o = rigid.oxygen_index;
        let h1 = rigid.hydrogen1_index;
        let h2 = rigid.hydrogen2_index;
        if o >= system.atoms.len() || h1 >= system.atoms.len() || h2 >= system.atoms.len() {
            return;
        }
        let m_o = system.atoms[o].mass;
        let m_h1 = system.atoms[h1].mass;
        let m_h2 = system.atoms[h2].mass;
        let total_mass = m_o + m_h1 + m_h2;
        if total_mass <= 0.0 {
            return;
        }

        let v_com = (system.atoms[o].velocity * m_o
            + system.atoms[h1].velocity * m_h1
            + system.atoms[h2].velocity * m_h2)
            / total_mass;
        let com = (system.atoms[o].position * m_o
            + system.atoms[h1].position * m_h1
            + system.atoms[h2].position * m_h2)
            / total_mass;
        let r_o = system.atoms[o].position - com;
        let r_h1 = system.atoms[h1].position - com;
        let r_h2 = system.atoms[h2].position - com;

        let l = m_o * r_o.cross(&(system.atoms[o].velocity - v_com))
            + m_h1 * r_h1.cross(&(system.atoms[h1].velocity - v_com))
            + m_h2 * r_h2.cross(&(system.atoms[h2].velocity - v_com));

        let inertia = |r: Vector3<f64>, m: f64| {
            let r2 = r.norm_squared();
            m * (Matrix3::identity() * r2 - r * r.transpose())
        };
        let i_tensor = inertia(r_o, m_o) + inertia(r_h1, m_h1) + inertia(r_h2, m_h2);
        let omega = i_tensor
            .try_inverse()
            .map(|inv| inv * l)
            .unwrap_or_else(Vector3::zeros);

        system.atoms[o].velocity = v_com + omega.cross(&r_o);
        system.atoms[h1].velocity = v_com + omega.cross(&r_h1);
        system.atoms[h2].velocity = v_com + omega.cross(&r_h2);
    }

    fn apply_constraint_positions(
        system: &mut System,
        constraints: Option<&[DistanceConstraint]>,
        options: &ConstraintOptions,
    ) {
        match options.mode {
            ConstraintMode::Flexible => {}
            ConstraintMode::ShakeRattle => {
                if let Some(c) = constraints {
                    shake_rattle::apply_shake(system, c, options.tolerance, options.max_iter);
                }
            }
            ConstraintMode::SettlePreferred => {
                if system.rigid_water.is_some() {
                    enforce_settle_position(system);
                } else if let Some(c) = constraints {
                    shake_rattle::apply_shake(system, c, options.tolerance, options.max_iter);
                }
            }
        }
    }

    fn apply_constraint_velocities(
        system: &mut System,
        constraints: Option<&[DistanceConstraint]>,
        options: &ConstraintOptions,
    ) {
        match options.mode {
            ConstraintMode::Flexible => {}
            ConstraintMode::ShakeRattle => {
                if let Some(c) = constraints {
                    shake_rattle::apply_rattle(system, c, options.tolerance, options.max_iter);
                }
            }
            ConstraintMode::SettlePreferred => {
                if system.rigid_water.is_some() {
                    enforce_settle_velocity(system);
                } else if let Some(c) = constraints {
                    shake_rattle::apply_rattle(system, c, options.tolerance, options.max_iter);
                }
            }
        }
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

            // --- half-step update for velocity and position ---
            // 1) velocity update (Verlet - half step)
            for (atom, a_o) in particles.iter_mut().zip(a_old.iter()) {
                atom.velocity += 0.5 * a_o * dt;
            }

            // 2) position update - needs to use the velocity that has been
            // updated using the half step method
            for atom in particles.iter_mut() {
                atom.update_position_verlet(dt); // the velocity has already been updated by a half step, so this will update the position by the half step
            }

            // ----

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
        run_md_nve_systems_with_constraints(
            systems,
            number_of_steps,
            dt,
            box_length,
            thermostat,
            None,
        );
    }

    pub fn run_md_nve_systems_with_constraints(
        systems: &mut Vec<System>,
        number_of_steps: i32,
        dt: f64,
        box_length: f64,
        thermostat: &str,
        constraint_options: Option<&ConstraintOptions>,
    ) {
        run_md_nve_systems_with_constraints_and_config(
            systems,
            number_of_steps,
            dt,
            box_length,
            thermostat,
            constraint_options,
            SystemSimulationConfig::default(),
        );
    }

    pub fn run_md_nve_systems_with_constraints_and_config(
        systems: &mut Vec<System>,
        number_of_steps: i32,
        dt: f64,
        box_length: f64,
        thermostat: &str,
        constraint_options: Option<&ConstraintOptions>,
        config: SystemSimulationConfig,
    ) {
        let mut values: Vec<f32> = Vec::new();
        let mut total_energy = 0.0;
        let mut kinetic_energy = 0.0;
        let mut potential_energy = 0.0;
        let mut nlist = SystemVerletNeighborList::new(
            config.cutoff,
            config.neighbor_skin,
            config.neighbor_rebuild_interval,
        );
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
        nlist.rebuild_pairs(systems, box_length);
        compute_intermolecular_forces_systems_with_pairs(
            systems,
            box_length,
            &nlist.pairs,
            config.cutoff,
        );
        let _ = add_electrostatic_forces_systems(systems, box_length, &config.pme);

        // this is only used if we apply nose hoover
        let mut xi_nose_hoover = vec![0.0; systems.len()];

        // --- time integration loop ---
        for _step in 0..number_of_steps {
            for (s, sys) in systems.iter_mut().enumerate() {
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
                if let Some(options) = constraint_options {
                    let constraints = options.constraints_by_system.get(s).map(|c| c.as_slice());
                    apply_constraint_positions(sys, constraints, options);
                }

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

            nlist.rebuild_if_needed(systems, box_length);
            compute_intermolecular_forces_systems_with_pairs(
                systems,
                box_length,
                &nlist.pairs,
                config.cutoff,
            );
            let electrostatic_energy =
                add_electrostatic_forces_systems(systems, box_length, &config.pme);

            for (s, sys) in systems.iter_mut().enumerate() {
                for a in sys.atoms.iter_mut() {
                    let a_new = a.force / a.mass;
                    a.update_velocity_verlet(a_new, dt);
                }
                if let Some(options) = constraint_options {
                    let constraints = options.constraints_by_system.get(s).map(|c| c.as_slice());
                    apply_constraint_velocities(sys, constraints, options);
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
            potential_energy +=
                intermolecular_site_site_energy_systems_cutoff(systems, box_length, config.cutoff);
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
                run_md_nve_particles(
                    particles,
                    number_of_steps,
                    dt,
                    box_length,
                    thermostat,
                    cutoff,
                );
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
    use nalgebra::Vector3;

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
        assert!(t.is_finite());
    }

    #[test]
    fn random_particle_initialization_avoids_close_contacts() {
        let state = lennard_jones_simulations::create_atoms_with_set_positions_and_velocities(
            25, 300.0, 30.0, 10.0, 10.0, false,
        )
        .expect("particle initialization should succeed");
        let lennard_jones_simulations::InitOutput::Particles(particles) = state else {
            panic!("expected particle initialization path");
        };

        let min_distance = particles
            .iter()
            .enumerate()
            .flat_map(|(i, a)| particles.iter().skip(i + 1).map(move |b| (a, b)))
            .map(|(a, b)| (a.position - b.position).norm())
            .fold(f64::INFINITY, f64::min);

        assert!(
            min_distance >= 0.9,
            "minimum sampled pair distance was {min_distance:.6}"
        );
    }

    #[test]
    fn exclusions_remove_intramolecular_electrostatics() {
        use crate::molecule::molecule::System;

        let mut sys = System::default();
        sys.atoms = vec![
            lennard_jones_simulations::Particle {
                id: 0,
                position: Vector3::new(0.0, 0.0, 0.0),
                velocity: Vector3::zeros(),
                force: Vector3::zeros(),
                lj_parameters: lennard_jones_simulations::LJParameters {
                    epsilon: 0.0,
                    sigma: 0.0,
                    number_of_atoms: 1,
                },
                mass: 1.0,
                energy: 0.0,
                atom_type: 0.0,
                charge: 1.0,
            },
            lennard_jones_simulations::Particle {
                id: 1,
                position: Vector3::new(1.0, 0.0, 0.0),
                velocity: Vector3::zeros(),
                force: Vector3::zeros(),
                lj_parameters: lennard_jones_simulations::LJParameters {
                    epsilon: 0.0,
                    sigma: 0.0,
                    number_of_atoms: 1,
                },
                mass: 1.0,
                energy: 0.0,
                atom_type: 0.0,
                charge: -1.0,
            },
        ];
        let mut systems = vec![sys.clone()];
        let e_before =
            lennard_jones_simulations::compute_electrostatic_forces_systems(&mut systems, 10.0);
        assert!(e_before.abs() > 1e-8);

        let mut excluded = sys;
        excluded.add_exclusion(0, 1);
        let mut systems_excluded = vec![excluded];
        let e_after = lennard_jones_simulations::compute_electrostatic_forces_systems(
            &mut systems_excluded,
            10.0,
        );
        let expected_self = -138.935_457_644 * 0.35 / std::f64::consts::PI.sqrt() * 2.0;
        assert!((e_after - expected_self).abs() < 1e-10);
        assert!(systems_excluded[0].atoms[0].force.norm() < 1e-12);
        assert!(systems_excluded[0].atoms[1].force.norm() < 1e-12);
    }

    fn test_particle(
        id: usize,
        position: Vector3<f64>,
        charge: f64,
    ) -> lennard_jones_simulations::Particle {
        lennard_jones_simulations::Particle {
            id,
            position,
            velocity: Vector3::zeros(),
            force: Vector3::zeros(),
            lj_parameters: lennard_jones_simulations::LJParameters {
                epsilon: 0.0,
                sigma: 0.0,
                number_of_atoms: 1,
            },
            mass: 1.0,
            energy: 0.0,
            atom_type: 0.0,
            charge,
        }
    }

    #[test]
    fn electrostatics_matches_analytic_two_charge_pair() {
        use crate::molecule::molecule::System;

        let mut sys = System::default();
        sys.atoms = vec![
            test_particle(0, Vector3::new(0.0, 0.0, 0.0), 1.0),
            test_particle(1, Vector3::new(1.0, 0.0, 0.0), -1.0),
        ];
        let mut systems = vec![sys];
        let no_screening = PmeConfig {
            alpha: 0.0,
            real_cutoff: 5.0,
            kmax: 0,
        };

        let energy = lennard_jones_simulations::compute_electrostatic_forces_systems_with_config(
            &mut systems,
            10.0,
            &no_screening,
        );
        let expected_energy = -138.935_457_644;
        assert!((energy - expected_energy).abs() < 1e-6);

        // Force on particle 0 should point toward +x with magnitude |k q1 q2 / r^2|.
        let expected_fx = 138.935_457_644;
        assert!((systems[0].atoms[0].force.x - expected_fx).abs() < 1e-6);
        assert!(systems[0].atoms[0].force.y.abs() < 1e-12);
        assert!(systems[0].atoms[0].force.z.abs() < 1e-12);
    }

    #[test]
    fn electrostatics_matches_direct_sum_for_small_tip3p_cluster() {
        use crate::molecule::molecule::System;

        // Two rigid TIP3P waters in nm with TIP3P charges (e).
        let mut sys = System::default();
        sys.atoms = vec![
            // water A
            test_particle(0, Vector3::new(0.00000, 0.00000, 0.00000), -0.834),
            test_particle(1, Vector3::new(0.09572, 0.00000, 0.00000), 0.417),
            test_particle(2, Vector3::new(-0.02399, 0.09266, 0.00000), 0.417),
            // water B
            test_particle(3, Vector3::new(0.30000, 0.10000, 0.05000), -0.834),
            test_particle(4, Vector3::new(0.39572, 0.10000, 0.05000), 0.417),
            test_particle(5, Vector3::new(0.27601, 0.19266, 0.05000), 0.417),
        ];
        let mut systems = vec![sys.clone()];
        let no_screening = PmeConfig {
            alpha: 0.0,
            real_cutoff: 9.0,
            kmax: 0,
        };
        let box_length = 4.0;

        let e_model = lennard_jones_simulations::compute_electrostatic_forces_systems_with_config(
            &mut systems,
            box_length,
            &no_screening,
        );

        let mut e_direct = 0.0;
        for i in 0..sys.atoms.len() {
            for j in (i + 1)..sys.atoms.len() {
                let dr = lennard_jones_simulations::minimum_image_convention(
                    sys.atoms[j].position - sys.atoms[i].position,
                    box_length,
                );
                let r = dr.norm();
                e_direct += 138.935_457_644 * sys.atoms[i].charge * sys.atoms[j].charge / r;
            }
        }

        let delta = (e_model - e_direct).abs();
        assert!(
            delta < 1e-3,
            "cluster electrostatic mismatch: model={e_model}, direct={e_direct}, delta={delta}"
        );
    }
}
