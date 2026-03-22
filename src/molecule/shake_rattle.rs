/*
SHAKE pseudocode:
-----------------

Given:
    positions r[i]
    velocities v[i]
    masses m[i]
    forces f[i]
    timestep dt
    constraints C = {(i, j, d_ij)}

For each MD step:

    # 1) Unconstrained integration step
    for each atom i:
        v_half[i] = v[i] + (dt/2) * f[i] / m[i]
        r_trial[i] = r[i] + dt * v_half[i]

    # 2) Iteratively enforce position constraints
    r_new = r_trial

    repeat until converged or max_iters reached:
        max_error = 0

        for each constraint (i, j, d):
            rij = r_new[i] - r_new[j]
            err = dot(rij, rij) - d*d

            max_error = max(max_error, abs(err))

            if abs(err) > tolerance:
                # solve for correction magnitude lambda
                # (details depend on integrator / formulation)
                lambda = approximate_constraint_multiplier(err, rij, m[i], m[j])

                # apply mass-weighted opposite corrections
                corr = lambda * rij
                r_new[i] = r_new[i] + corr / m[i]
                r_new[j] = r_new[j] - corr / m[j]

        if max_error < tolerance:
            break

    # 3) Recompute forces from constrained positions
    f_new = forces(r_new)

    # 4) Finish velocity update
    for each atom i:
        v_new[i] = v_half[i] + (dt/2) * f_new[i] / m[i]



Rattle pseudocode
-----------------

RATTLE does the same position correction, but then adds a second stage so that the velocities are also consistent with the constraints.

For a fixed bond distance, not only must

∣
𝑟
𝑖
−
𝑟
𝑗
∣
2
=
𝑑
2
∣r
i
​
 −r
j
​
 ∣
2
 =d
2

hold, but also its time derivative must be zero, which means there should be no relative velocity along the bond direction. LAMMPS explicitly notes that fix rattle modifies forces and velocities, while OpenMM separately provides “Constrain Velocities.”


Given:
    positions r[i]
    velocities v[i]
    masses m[i]
    forces f[i]
    timestep dt
    constraints C = {(i, j, d_ij)}

For each MD step:

    # 1) Unconstrained integration step
    for each atom i:
        v_half[i] = v[i] + (dt/2) * f[i] / m[i]
        r_trial[i] = r[i] + dt * v_half[i]

    # 2) Iteratively enforce position constraints
    r_new = r_trial

    repeat until converged or max_iters reached:
        max_error = 0

        for each constraint (i, j, d):
            rij = r_new[i] - r_new[j]
            err = dot(rij, rij) - d*d

            max_error = max(max_error, abs(err))

            if abs(err) > tolerance:
                # solve for correction magnitude lambda
                # (details depend on integrator / formulation)
                lambda = approximate_constraint_multiplier(err, rij, m[i], m[j])

                # apply mass-weighted opposite corrections
                corr = lambda * rij
                r_new[i] = r_new[i] + corr / m[i]
                r_new[j] = r_new[j] - corr / m[j]

        if max_error < tolerance:
            break

    # 3) Recompute forces from constrained positions
    f_new = forces(r_new)

    # 4) Finish velocity update
    for each atom i:
        v_new[i] = v_half[i] + (dt/2) * f_new[i] / m[i]

*/
