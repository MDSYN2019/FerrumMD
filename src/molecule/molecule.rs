/*

Expanding your rust based molecular dynamics (MD) simulation from point particles
to molecules with bonded interactions and force fields requires adding new types
of interactions

 */

use crate::lennard_jones_simulations::minimum_image_convention;
use crate::lennard_jones_simulations::InitOutput;
use crate::lennard_jones_simulations::LJParameters;
use crate::lennard_jones_simulations::Particle;

use nalgebra::Vector3;
use std::collections::{HashMap, HashSet};

#[derive(Copy, Clone)]
pub struct SimpleBond {
    pub i: usize,
    pub j: usize,
    pub k: f64,
    pub r0: f64,
}

#[derive(Clone)]
pub struct Atom {
    // Definition of the atom, force on the atom, as well as its current position and velocity
    // as well as the mass and charge
    pub id: usize,
    pub position: Vector3<f64>,
    pub velocity: Vector3<f64>,
    pub force: Vector3<f64>,
    pub atom_type: usize,
    pub mass: f64,
    pub charge: f64,
}

#[derive(Clone, Debug)]
pub struct Bond {
    pub atom1: usize,
    pub atom2: usize,
    pub k: f64,
    pub r0: f64,
}

#[derive(Clone, Debug)]
pub struct Angle {
    pub atom1: usize,
    pub atom2: usize,
    pub atom3: usize,
    pub k: f64,
    pub theta0: f64,
}

#[derive(Clone, Debug)]
pub struct Dihedral {
    pub atom1: usize,
    pub atom2: usize,
    pub atom3: usize,
    pub atom4: usize,
    pub k: f64,
    pub multiplicity: usize,
    pub phase: f64,
}

#[derive(Clone, Debug)]
pub struct Improper {
    pub atom1: usize,
    pub atom2: usize,
    pub atom3: usize,
    pub atom4: usize,
    pub k: f64,
    pub psi0: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct NonbondedPairScaling {
    pub lj_scale: f64,
    pub coulomb_scale: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct RigidWaterDefinition {
    pub oxygen_index: usize,
    pub hydrogen1_index: usize,
    pub hydrogen2_index: usize,
    pub oh_distance: f64,
    pub hh_distance: f64,
}

#[derive(Copy, Clone)]
pub struct NonBondedType {
    pub mass: f64,
    pub charge: f64,
    pub sigma: f64,
    pub epsilon: f64,
}

pub struct ForceField {
    //pub atom_types: HashMap<String, AtomTypeParams>,
    pub bond_types: Vec<Bond>,
    pub angle_types: Vec<Angle>,
}

pub fn bond_distance(a1: &Atom, a2: &Atom) -> f64 {
    // Compute the bond distance between atoms
    (a1.position - a2.position).norm()
}

#[derive(Default, Clone)]
pub struct MoleculeTemplate {
    pub name: String,
    pub atom_types: Vec<String>,      // len
    pub positions: Vec<Vector3<f64>>, // x y z for each atom
    pub bonds: Vec<(usize, usize, f64, f64)>,
    pub exclusion_1_4_scale: Option<f64>, // (i, j, k, k_theta, theta_0)
}

#[derive(Clone, Default, Debug)]
pub struct System {
    pub atoms: Vec<Particle>,
    pub bonds: Vec<Bond>,
    pub angles: Vec<Angle>,
    pub dihedrals: Vec<Dihedral>,
    pub impropers: Vec<Improper>,
    pub exclusions: HashSet<(usize, usize)>,
    pub pair_scalings: HashMap<(usize, usize), NonbondedPairScaling>,
    pub rigid_water: Option<RigidWaterDefinition>,
}

// System is all the atoms (global), bonded terms in global indices, and exclusion sets
impl System {
    fn canonical_pair(i: usize, j: usize) -> (usize, usize) {
        if i <= j {
            (i, j)
        } else {
            (j, i)
        }
    }

    pub fn add_exclusion(&mut self, i: usize, j: usize) {
        self.exclusions.insert(Self::canonical_pair(i, j));
    }

    pub fn set_pair_scaling(&mut self, i: usize, j: usize, scaling: NonbondedPairScaling) {
        self.pair_scalings
            .insert(Self::canonical_pair(i, j), scaling);
    }

    pub fn nonbonded_scaling_for_pair(&self, i: usize, j: usize) -> Option<NonbondedPairScaling> {
        if i == j || self.exclusions.contains(&Self::canonical_pair(i, j)) {
            return None;
        }
        self.pair_scalings
            .get(&Self::canonical_pair(i, j))
            .copied()
            .or(Some(NonbondedPairScaling {
                lj_scale: 1.0,
                coulomb_scale: 1.0,
            }))
    }
}

pub fn compute_bond_force(atoms: &mut Vec<Particle>, bond: &Bond, box_length: f64) -> f64 {
    /*
    Compute the bond energy,
     */
    let (i, j) = (bond.atom1, bond.atom2); // get atoms#
    let r_vec = atoms[j].position - atoms[i].position; // get vector for position
    let rij_mic = minimum_image_convention(r_vec, box_length);
    let r = rij_mic.norm(); // get distance
    let dr = r - bond.r0; // the difference between the current position and the equilibrium position
    let f_mag = -bond.k * dr; // force magnitude
    let f_vec = (r_vec / r) * f_mag;

    atoms[i].force += f_vec;
    atoms[j].force -= f_vec;

    0.5 * bond.k * dr * dr // return the bond energy
}

pub fn compute_electostatic_bond_short_force(atoms: &mut Vec<Particle>, _box_length: f64) -> f64 {
    /*
    Compute the short range real space component of the electrostatic interaction

    https://computecanada.github.io/molmodsim-md-theory-lesson-novice/06-electrostatics/index.html - useful link

    Computing Coulomb potenials is often the most time consuming part of any MD simulation

     */
    let mut total_short_range_potential = 0.0;
    let e_0 = 1.0;
    for i in 0..atoms.len() {
        for j in (i + 1)..atoms.len() {
            // This needs to be properly represent the coloumbing potential - this is a crappy dummy at the moment
            total_short_range_potential += ((atoms[i].charge * atoms[j].charge)
                / (4.0 * 3.14 * e_0))
                / (atoms[0].position - atoms[1].position).norm()
        }
    }
    total_short_range_potential
}

fn angle_value(atoms: &[Particle], angle: &Angle, box_length: f64) -> f64 {
    let r21 = minimum_image_convention(
        atoms[angle.atom1].position - atoms[angle.atom2].position,
        box_length,
    );
    let r23 = minimum_image_convention(
        atoms[angle.atom3].position - atoms[angle.atom2].position,
        box_length,
    );

    let n1 = r21.norm();
    let n2 = r23.norm();
    if n1 <= 1e-12 || n2 <= 1e-12 {
        return angle.theta0;
    }

    let cos_theta = (r21.dot(&r23) / (n1 * n2)).clamp(-1.0, 1.0);
    cos_theta.acos()
}

fn dihedral_value(atoms: &[Particle], dihedral: &Dihedral, box_length: f64) -> f64 {
    let b1 = minimum_image_convention(
        atoms[dihedral.atom2].position - atoms[dihedral.atom1].position,
        box_length,
    );
    let b2 = minimum_image_convention(
        atoms[dihedral.atom3].position - atoms[dihedral.atom2].position,
        box_length,
    );
    let b3 = minimum_image_convention(
        atoms[dihedral.atom4].position - atoms[dihedral.atom3].position,
        box_length,
    );

    let n1 = b1.cross(&b2);
    let n2 = b2.cross(&b3);
    let b2_norm = b2.norm();
    if n1.norm() <= 1e-12 || n2.norm() <= 1e-12 || b2_norm <= 1e-12 {
        return 0.0;
    }

    let b2_hat = b2 / b2_norm;
    let m1 = n1.cross(&b2_hat);

    let x = n1.dot(&n2);
    let y = m1.dot(&n2);
    y.atan2(x)
}

fn improper_value(atoms: &[Particle], improper: &Improper, box_length: f64) -> f64 {
    let as_dihedral = Dihedral {
        atom1: improper.atom1,
        atom2: improper.atom2,
        atom3: improper.atom3,
        atom4: improper.atom4,
        k: improper.k,
        multiplicity: 1,
        phase: 0.0,
    };
    dihedral_value(atoms, &as_dihedral, box_length)
}

pub fn compute_angle_force(atoms: &mut [Particle], angle: &Angle, box_length: f64) -> f64 {
    let theta = angle_value(atoms, angle, box_length);
    let dtheta = theta - angle.theta0;
    let energy = 0.5 * angle.k * dtheta * dtheta;

    let atom_indices = [angle.atom1, angle.atom2, angle.atom3];
    let h = 1e-6;

    for &idx in &atom_indices {
        for dim in 0..3 {
            atoms[idx].position[dim] += h;
            let e_plus =
                0.5 * angle.k * (angle_value(atoms, angle, box_length) - angle.theta0).powi(2);
            atoms[idx].position[dim] -= 2.0 * h;
            let e_minus =
                0.5 * angle.k * (angle_value(atoms, angle, box_length) - angle.theta0).powi(2);
            atoms[idx].position[dim] += h;

            let d_e = (e_plus - e_minus) / (2.0 * h);
            atoms[idx].force[dim] += -d_e;
        }
    }

    energy
}

pub fn compute_dihedral_force(atoms: &mut [Particle], dihedral: &Dihedral, box_length: f64) -> f64 {
    let phi = dihedral_value(atoms, dihedral, box_length);
    let n = dihedral.multiplicity as f64;
    let energy = dihedral.k * (1.0 + (n * phi - dihedral.phase).cos());

    let atom_indices = [
        dihedral.atom1,
        dihedral.atom2,
        dihedral.atom3,
        dihedral.atom4,
    ];
    let h = 1e-6;

    for &idx in &atom_indices {
        for dim in 0..3 {
            atoms[idx].position[dim] += h;
            let e_plus = dihedral.k
                * (1.0
                    + ((n * dihedral_value(atoms, dihedral, box_length)) - dihedral.phase).cos());
            atoms[idx].position[dim] -= 2.0 * h;
            let e_minus = dihedral.k
                * (1.0
                    + ((n * dihedral_value(atoms, dihedral, box_length)) - dihedral.phase).cos());
            atoms[idx].position[dim] += h;

            let d_e = (e_plus - e_minus) / (2.0 * h);
            atoms[idx].force[dim] += -d_e;
        }
    }

    energy
}

pub fn compute_improper_force(atoms: &mut [Particle], improper: &Improper, box_length: f64) -> f64 {
    let psi = improper_value(atoms, improper, box_length);
    let dpsi = psi - improper.psi0;
    let energy = 0.5 * improper.k * dpsi * dpsi;

    let atom_indices = [
        improper.atom1,
        improper.atom2,
        improper.atom3,
        improper.atom4,
    ];
    let h = 1e-6;

    for &idx in &atom_indices {
        for dim in 0..3 {
            atoms[idx].position[dim] += h;
            let e_plus = 0.5
                * improper.k
                * (improper_value(atoms, improper, box_length) - improper.psi0).powi(2);
            atoms[idx].position[dim] -= 2.0 * h;
            let e_minus = 0.5
                * improper.k
                * (improper_value(atoms, improper, box_length) - improper.psi0).powi(2);
            atoms[idx].position[dim] += h;

            let d_e = (e_plus - e_minus) / (2.0 * h);
            atoms[idx].force[dim] += -d_e;
        }
    }

    energy
}

pub fn apply_all_bonded_forces_and_energy(
    atoms: &mut Vec<Particle>,
    bonds: &[Bond],
    angles: &[Angle],
    dihedrals: &[Dihedral],
    impropers: &[Improper],
    box_length: f64,
) -> f64 {
    let mut energy = 0.0;

    for b in bonds {
        energy += compute_bond_force(atoms, b, box_length);
    }
    for angle in angles {
        energy += compute_angle_force(atoms, angle, box_length);
    }
    for dihedral in dihedrals {
        energy += compute_dihedral_force(atoms, dihedral, box_length);
    }
    for improper in impropers {
        energy += compute_improper_force(atoms, improper, box_length);
    }

    energy
}

pub fn apply_bonded_forces_and_energy(
    atoms: &mut Vec<Particle>,
    bonds: &[Bond],
    box_length: f64,
) -> f64 {
    apply_all_bonded_forces_and_energy(atoms, bonds, &[], &[], &[], box_length)
}

pub fn make_h2_system() -> System {
    /*
    Reduced units:

    mass = 1.0 for each H (you can sue 1.0 amu reduced)


    The H-H distance r oscillates around r0 = 0.74

     */
    let r0 = 0.74;
    let k = 100.0;
    let x = 0.5 * r0;

    let mut atoms = vec![
        Particle {
            id: 0,
            position: Vector3::new(-x, 0.0, 0.0),
            velocity: Vector3::new(1.0, 1.0, 1.0),
            force: Vector3::zeros(),
            atom_type: 0.0,
            mass: 1.0,
            charge: 0.1,
            energy: 0.0,
            lj_parameters: (LJParameters {
                epsilon: 1.0,
                sigma: 1.0,
                number_of_atoms: 3, // this needs to be corrected
            }),
        },
        Particle {
            id: 1,
            position: Vector3::new(x, 0.0, 0.0),
            velocity: Vector3::new(1.0, -1.0, 1.0),
            force: Vector3::zeros(),
            atom_type: 0.0,
            mass: 1.0,
            charge: -0.1,
            energy: 0.0,
            lj_parameters: (LJParameters {
                epsilon: 1.0,
                sigma: 1.0,
                number_of_atoms: 3, // this needs to be corrected
            }),
        },
    ];

    let stretch = 0.05;
    atoms[1].position.x += 0.5 * stretch;
    atoms[0].position.x -= 0.5 * stretch;

    let bonds = vec![Bond {
        atom1: 0,
        atom2: 1,
        k,
        r0,
    }];

    System {
        atoms,
        bonds,
        angles: vec![],
        dihedrals: vec![],
        impropers: vec![],
        ..System::default()
    }
}

pub fn create_systems(system: &System, number_of_molecules: i32) -> InitOutput {
    /*
    Create n number of particles
     */
    let mut molecules: Vec<System> = Vec::new();

    if number_of_molecules <= 0 {
        return InitOutput::Systems(molecules);
    }

    // Build a simple cubic placement to avoid placing every molecule at the
    // exact same coordinates (which causes extreme LJ overlaps/energies).
    let mut min = Vector3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut max = Vector3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for atom in &system.atoms {
        min.x = min.x.min(atom.position.x);
        min.y = min.y.min(atom.position.y);
        min.z = min.z.min(atom.position.z);
        max.x = max.x.max(atom.position.x);
        max.y = max.y.max(atom.position.y);
        max.z = max.z.max(atom.position.z);
    }

    let span = max - min;
    let base_spacing = span.x.max(span.y).max(span.z);
    let spacing = (base_spacing + 2.5).max(3.0);
    let n_side = (number_of_molecules as f64).cbrt().ceil() as i32;

    for m in 0..number_of_molecules {
        let ix = m % n_side;
        let iy = (m / n_side) % n_side;
        let iz = m / (n_side * n_side);

        let mut sys_clone = system.clone();
        let shift = Vector3::new(
            ix as f64 * spacing,
            iy as f64 * spacing,
            iz as f64 * spacing,
        );

        for atom in &mut sys_clone.atoms {
            atom.position += shift;
        }

        molecules.push(sys_clone);
    }

    // output as the enum we want which will be a valid input to run_md_nve
    InitOutput::Systems(molecules)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    #[test]
    fn test_atom() {
        let test_atom = Atom {
            id: 1,
            position: Vector3::new(1.0, 2.0, 3.0),
            velocity: Vector3::new(0.0, 0.0, 0.0),
            force: Vector3::new(0.0, 0.0, 0.0),
            atom_type: 0,
            mass: 12.01,
            charge: -3.0,
        };
        assert_eq!(test_atom.id, 1);
    }

    #[test]
    fn test_create_systems_offsets_molecules() {
        let h2 = make_h2_system();
        let created = create_systems(&h2, 8);

        let systems = match created {
            InitOutput::Systems(systems) => systems,
            InitOutput::Particles(_) => panic!("expected systems output"),
        };

        assert_eq!(systems.len(), 8);

        let first = systems[0].atoms[0].position;
        let mut distinct_positions = 0;
        for sys in systems.iter().skip(1) {
            if (sys.atoms[0].position - first).norm() > 1e-12 {
                distinct_positions += 1;
            }
        }

        assert!(
            distinct_positions > 0,
            "cloned molecules should not all occupy the same coordinates"
        );
    }

    #[test]
    fn test_bond_distance_happy() {
        let atom1 = Atom {
            id: 0,
            position: Vector3::new(0.0, 0.0, 0.0),
            velocity: Vector3::zeros(),
            force: Vector3::zeros(),
            atom_type: 1,
            mass: 1.0,
            charge: 0.0,
        };

        let atom2 = Atom {
            id: 1,
            position: Vector3::new(3.0, 0.0, 0.0),
            ..atom1.clone()
        };

        let dist = bond_distance(&atom1, &atom2);
        assert!((dist - 3.0).abs() < 1e-6) // happy path
    }
    #[test]
    fn test_angle_force_energy_zero_at_equilibrium() {
        let mut atoms = vec![
            Particle {
                id: 0,
                position: Vector3::new(1.0, 0.0, 0.0),
                velocity: Vector3::zeros(),
                force: Vector3::zeros(),
                atom_type: 0.0,
                mass: 1.0,
                charge: 0.0,
                energy: 0.0,
                lj_parameters: LJParameters {
                    epsilon: 1.0,
                    sigma: 1.0,
                    number_of_atoms: 1,
                },
            },
            Particle {
                id: 1,
                position: Vector3::new(0.0, 0.0, 0.0),
                velocity: Vector3::zeros(),
                force: Vector3::zeros(),
                atom_type: 0.0,
                mass: 1.0,
                charge: 0.0,
                energy: 0.0,
                lj_parameters: LJParameters {
                    epsilon: 1.0,
                    sigma: 1.0,
                    number_of_atoms: 1,
                },
            },
            Particle {
                id: 2,
                position: Vector3::new(0.0, 1.0, 0.0),
                velocity: Vector3::zeros(),
                force: Vector3::zeros(),
                atom_type: 0.0,
                mass: 1.0,
                charge: 0.0,
                energy: 0.0,
                lj_parameters: LJParameters {
                    epsilon: 1.0,
                    sigma: 1.0,
                    number_of_atoms: 1,
                },
            },
        ];

        let angle = Angle {
            atom1: 0,
            atom2: 1,
            atom3: 2,
            k: 10.0,
            theta0: std::f64::consts::FRAC_PI_2,
        };

        let e = compute_angle_force(&mut atoms, &angle, 10.0);
        assert!(e.abs() < 1e-8);
    }

    #[test]
    fn test_dihedral_energy_phase_shift() {
        let atoms = vec![
            Particle {
                id: 0,
                position: Vector3::new(0.0, 0.0, 0.0),
                velocity: Vector3::zeros(),
                force: Vector3::zeros(),
                atom_type: 0.0,
                mass: 1.0,
                charge: 0.0,
                energy: 0.0,
                lj_parameters: LJParameters {
                    epsilon: 1.0,
                    sigma: 1.0,
                    number_of_atoms: 1,
                },
            },
            Particle {
                id: 1,
                position: Vector3::new(1.0, 0.0, 0.0),
                velocity: Vector3::zeros(),
                force: Vector3::zeros(),
                atom_type: 0.0,
                mass: 1.0,
                charge: 0.0,
                energy: 0.0,
                lj_parameters: LJParameters {
                    epsilon: 1.0,
                    sigma: 1.0,
                    number_of_atoms: 1,
                },
            },
            Particle {
                id: 2,
                position: Vector3::new(1.0, 1.0, 0.0),
                velocity: Vector3::zeros(),
                force: Vector3::zeros(),
                atom_type: 0.0,
                mass: 1.0,
                charge: 0.0,
                energy: 0.0,
                lj_parameters: LJParameters {
                    epsilon: 1.0,
                    sigma: 1.0,
                    number_of_atoms: 1,
                },
            },
            Particle {
                id: 3,
                position: Vector3::new(1.0, 1.0, 1.0),
                velocity: Vector3::zeros(),
                force: Vector3::zeros(),
                atom_type: 0.0,
                mass: 1.0,
                charge: 0.0,
                energy: 0.0,
                lj_parameters: LJParameters {
                    epsilon: 1.0,
                    sigma: 1.0,
                    number_of_atoms: 1,
                },
            },
        ];

        let dih = Dihedral {
            atom1: 0,
            atom2: 1,
            atom3: 2,
            atom4: 3,
            k: 2.0,
            multiplicity: 1,
            phase: 0.0,
        };

        let phi = dihedral_value(&atoms, &dih, 10.0);
        let e = dih.k * (1.0 + (phi - dih.phase).cos());
        assert!(e.is_finite());
    }
}
