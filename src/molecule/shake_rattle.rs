pub mod shake_rattle {
    use crate::lennard_jones_simulations::Particle;
    use crate::molecule::molecule::System;

    fn rattle(
        tolerance: f64,
        system: &mut System,
        index_i: usize,
        index_j: usize,
        target_distance: f64,
        max_iter: usize,
    ) -> () {
        /*
        Rattle
        ------

        Rattle is basicall a velocity/verlet + constraint correction

        -> first correct positions so bond lengths are satisifed
        -> then correct the velocities so the velocity constraints are satisfied

        Essentially, we are nudging the atoms along the bond direction until the distance matches the target

        so:

        -> if the bond length is too long, pull them together
        -> if the bond is too short, push them apart

         */

        // compute the inverse masses
        let inv_mi = 1.0 / system.atoms[index_i].mass;
        let inv_mj = 1.0 / system.atoms[index_j].mass;

        //let mut lambda_ij = 0.0; //
        for _ in 0..max_iter {
            let r_vec = system.atoms[index_j].position - system.atoms[index_i].position; // get vector for position
            let dist_sq = r_vec.dot(&r_vec); // get the distance squared

            // position correction - for each constrained bond (i,j) with
            // target distance (target_distance)

            let mut constraint_error = dist_sq - target_distance.powi(2);
            // we must iteratively adjust positions until constraint_error is roughly 0

            if constraint_error.abs() > tolerance {
                break;
            }

            // avoid divide-by-zero
            if dist_sq < 1e-12 {
                panic!("Bond distance is too small during RATTLE correction");
            }

            // approximate the scalar correction
            let lambda_ij = -constraint_error / (2.0 * (inv_mi + inv_mj) * dist_sq);

            let correction = lambda_ij * r_vec;

            // Apply mass-weighted corrections
            system.atoms[index_i].position -= correction * inv_mi;
            system.atoms[index_j].position += correction * inv_mj;
        }
    }

    fn shake(tolerance: f64) -> () {}
}
