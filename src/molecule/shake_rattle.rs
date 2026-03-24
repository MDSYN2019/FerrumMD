/*
SHAKE: fixes positions so constraints like bond lengths are satisfied

RATTLE: does SHAKE plus an extra correction for velocities

*/

pub mod shake_rattle {
    use nalgebra::Vector3;

    use crate::molecule::molecule::System;

    pub fn tip3p_constraints_from_system(system: &System) -> Result<Vec<(usize, usize, f64)>, String> {
        if system.atoms.len() != 3 {
            return Err("TIP3P constraint builder expects exactly 3 atoms per molecule.".to_string());
        }

        let oh1 = (system.atoms[1].position - system.atoms[0].position).norm();
        let oh2 = (system.atoms[2].position - system.atoms[0].position).norm();
        let hh = (system.atoms[2].position - system.atoms[1].position).norm();

        Ok(vec![(0, 1, oh1), (0, 2, oh2), (1, 2, hh)])
    }

    fn shake_bond(
        tolerance: f64,
        system: &mut System,
        index_i: usize,
        index_j: usize,
        target_distance: f64,
        max_iter: usize,
    ) {
        let inv_mi = 1.0 / system.atoms[index_i].mass;
        let inv_mj = 1.0 / system.atoms[index_j].mass;

        for _ in 0..max_iter {
            let r_vec = system.atoms[index_j].position - system.atoms[index_i].position;
            let dist_sq = r_vec.dot(&r_vec);
            let constraint_error = dist_sq - target_distance.powi(2);

            if constraint_error.abs() < tolerance {
                break;
            }

            if dist_sq <= 1e-12 {
                break;
            }
            let denominator = 2.0 * (inv_mi + inv_mj) * dist_sq;

            if denominator <= 1e-12 {
                break;
            }

            let lambda_ij = -constraint_error / denominator;
            let correction = lambda_ij * r_vec;

            system.atoms[index_i].position -= correction * inv_mi;
            system.atoms[index_j].position += correction * inv_mj;
        }
    }

    fn rattle_bond(
        tolerance: f64,
        system: &mut System,
        index_i: usize,
        index_j: usize,
        max_iter: usize,
    ) {
        let inv_mi = 1.0 / system.atoms[index_i].mass;
        let inv_mj = 1.0 / system.atoms[index_j].mass;

        for _ in 0..max_iter {
            let r_vec: Vector3<f64> = system.atoms[index_j].position - system.atoms[index_i].position;
            let v_vec: Vector3<f64> = system.atoms[index_j].velocity - system.atoms[index_i].velocity;
            let dist_sq = r_vec.dot(&r_vec);
            if dist_sq <= 1e-12 {
                break;
            }

            let constraint_velocity_error = r_vec.dot(&v_vec);
            if constraint_velocity_error.abs() < tolerance {
                break;
            }

            let denominator = (inv_mi + inv_mj) * dist_sq;
            if denominator <= 1e-12 {
                break;
            }

            let mu = -constraint_velocity_error / denominator;
            let velocity_correction = mu * r_vec;

            system.atoms[index_i].velocity -= velocity_correction * inv_mi;
            system.atoms[index_j].velocity += velocity_correction * inv_mj;
        }
    }

    pub fn apply_shake(
        system: &mut System,
        constraints: &[(usize, usize, f64)],
        tolerance: f64,
        max_iter: usize,
    ) {
        for _ in 0..max_iter {
            let mut converged = true;
            for &(i, j, d_ij) in constraints {
                let r_vec = system.atoms[j].position - system.atoms[i].position;
                let error = r_vec.dot(&r_vec) - d_ij.powi(2);
                if error.abs() >= tolerance {
                    converged = false;
                }
                shake_bond(tolerance, system, i, j, d_ij, max_iter);
            }
            if converged {
                break;
            }
        }
    }

    pub fn apply_rattle(
        system: &mut System,
        constraints: &[(usize, usize, f64)],
        tolerance: f64,
        max_iter: usize,
    ) {
        for _ in 0..max_iter {
            let mut converged = true;
            for &(i, j, _d_ij) in constraints {
                let r_vec = system.atoms[j].position - system.atoms[i].position;
                let v_vec = system.atoms[j].velocity - system.atoms[i].velocity;
                if r_vec.dot(&v_vec).abs() >= tolerance {
                    converged = false;
                }
                rattle_bond(tolerance, system, i, j, max_iter);
            }
            if converged {
                break;
            }
        }
    }
}
