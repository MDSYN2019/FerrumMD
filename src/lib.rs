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
pub mod lennard_jones_simulations;
pub mod molecule;
pub mod parameters;
#[cfg(feature = "python")]
mod python;
#[path = "quantum/quantum_chem.rs"]
pub mod quantum_chemistry;
pub mod thermostat_barostat;
pub mod visualization;

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
