use crate::lennard_jones_simulations::Particle;
use nalgebra::Vector3;
use std::sync::Arc;
use std::thread;

#[inline]
fn lj_force_prefactor_from_r2(r2: f64, sigma: f64, epsilon: f64) -> f64 {
    if r2 <= 1e-24 {
        return 0.0;
    }
    let inv_r2 = 1.0 / r2;
    let sigma2 = sigma * sigma;
    let sr2 = sigma2 * inv_r2;
    let sr6 = sr2 * sr2 * sr2;
    let sr12 = sr6 * sr6;
    24.0 * epsilon * (2.0 * sr12 - sr6) * inv_r2
}

pub fn compute_forces_python_baseline(particles: &[Particle]) -> Vec<Vector3<f64>> {
    let n = particles.len();
    let mut forces = vec![Vector3::zeros(); n];

    for i in 0..n {
        let mut fi = Vector3::zeros();
        for j in 0..n {
            if i == j {
                continue;
            }
            let dr = particles[j].position - particles[i].position;
            let r2 = dr.dot(&dr);
            if r2 <= 1e-24 {
                continue;
            }
            let r = r2.sqrt();
            let sigma = 0.5 * (particles[i].lj_parameters.sigma + particles[j].lj_parameters.sigma);
            let epsilon =
                (particles[i].lj_parameters.epsilon * particles[j].lj_parameters.epsilon).sqrt();

            // Intentionally Python-like scalar formulation with powf-heavy arithmetic.
            let sr6 = (sigma / r).powf(6.0);
            let sr12 = sr6 * sr6;
            let f_mag = 24.0 * epsilon * (2.0 * sr12 - sr6) / r;
            fi -= (dr / r) * f_mag;
        }
        forces[i] = fi;
    }

    forces
}

pub fn compute_forces_simd_parallel(particles: &[Particle]) -> Vec<Vector3<f64>> {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if std::is_x86_feature_detected!("avx") {
            return compute_forces_simd_parallel_avx(particles);
        }
    }

    compute_forces_parallel_scalar(particles)
}

fn compute_forces_parallel_scalar(particles: &[Particle]) -> Vec<Vector3<f64>> {
    let n = particles.len();
    let threads = thread::available_parallelism().map_or(1, |n| n.get());
    let chunk = n.div_ceil(threads);
    let particles = Arc::new(particles.to_vec());
    let mut handles = Vec::new();

    for t in 0..threads {
        let start = t * chunk;
        if start >= n {
            break;
        }
        let end = ((t + 1) * chunk).min(n);
        let particles_ref = Arc::clone(&particles);

        handles.push(thread::spawn(move || {
            let mut out = vec![Vector3::zeros(); end - start];
            for (local_i, i) in (start..end).enumerate() {
                let pi = &particles_ref[i];
                let mut fi = Vector3::zeros();
                for (j, pj) in particles_ref.iter().enumerate() {
                    if i == j {
                        continue;
                    }
                    let dr = pj.position - pi.position;
                    let r2 = dr.dot(&dr);
                    let sigma = 0.5 * (pi.lj_parameters.sigma + pj.lj_parameters.sigma);
                    let epsilon = (pi.lj_parameters.epsilon * pj.lj_parameters.epsilon).sqrt();
                    let pref = lj_force_prefactor_from_r2(r2, sigma, epsilon);
                    fi -= dr * pref;
                }
                out[local_i] = fi;
            }
            (start, out)
        }));
    }

    let mut forces = vec![Vector3::zeros(); n];
    for handle in handles {
        let (start, block) = handle.join().expect("force worker thread panicked");
        for (k, f) in block.into_iter().enumerate() {
            forces[start + k] = f;
        }
    }

    forces
}

#[cfg(target_arch = "x86_64")]
fn compute_forces_simd_parallel_avx(particles: &[Particle]) -> Vec<Vector3<f64>> {
    use std::arch::x86_64::*;

    let n = particles.len();
    let xs: Vec<f64> = particles.iter().map(|p| p.position.x).collect();
    let ys: Vec<f64> = particles.iter().map(|p| p.position.y).collect();
    let zs: Vec<f64> = particles.iter().map(|p| p.position.z).collect();
    let sigmas: Vec<f64> = particles.iter().map(|p| p.lj_parameters.sigma).collect();
    let eps: Vec<f64> = particles.iter().map(|p| p.lj_parameters.epsilon).collect();

    let threads = thread::available_parallelism().map_or(1, |n| n.get());
    let chunk = n.div_ceil(threads);

    let xs = Arc::new(xs);
    let ys = Arc::new(ys);
    let zs = Arc::new(zs);
    let sigmas = Arc::new(sigmas);
    let eps = Arc::new(eps);

    let mut handles = Vec::new();

    for t in 0..threads {
        let start = t * chunk;
        if start >= n {
            break;
        }
        let end = ((t + 1) * chunk).min(n);

        let xs_ref = Arc::clone(&xs);
        let ys_ref = Arc::clone(&ys);
        let zs_ref = Arc::clone(&zs);
        let sigmas_ref = Arc::clone(&sigmas);
        let eps_ref = Arc::clone(&eps);

        handles.push(thread::spawn(move || {
            let mut out = vec![Vector3::zeros(); end - start];

            for (local_i, i) in (start..end).enumerate() {
                let mut fx = 0.0;
                let mut fy = 0.0;
                let mut fz = 0.0;

                let xi = xs_ref[i];
                let yi = ys_ref[i];
                let zi = zs_ref[i];
                let si = sigmas_ref[i];
                let ei = eps_ref[i];

                let mut j = 0usize;
                while j + 4 <= n {
                    unsafe {
                        let xj = _mm256_loadu_pd(xs_ref[j..].as_ptr());
                        let yj = _mm256_loadu_pd(ys_ref[j..].as_ptr());
                        let zj = _mm256_loadu_pd(zs_ref[j..].as_ptr());
                        let sj = _mm256_loadu_pd(sigmas_ref[j..].as_ptr());
                        let ej = _mm256_loadu_pd(eps_ref[j..].as_ptr());

                        let xi_v = _mm256_set1_pd(xi);
                        let yi_v = _mm256_set1_pd(yi);
                        let zi_v = _mm256_set1_pd(zi);

                        let dx = _mm256_sub_pd(xj, xi_v);
                        let dy = _mm256_sub_pd(yj, yi_v);
                        let dz = _mm256_sub_pd(zj, zi_v);

                        let r2 = _mm256_add_pd(
                            _mm256_mul_pd(dx, dx),
                            _mm256_add_pd(_mm256_mul_pd(dy, dy), _mm256_mul_pd(dz, dz)),
                        );

                        let inv_r2 = _mm256_div_pd(_mm256_set1_pd(1.0), r2);
                        let sigma = _mm256_mul_pd(
                            _mm256_set1_pd(0.5),
                            _mm256_add_pd(_mm256_set1_pd(si), sj),
                        );
                        let epsilon = _mm256_sqrt_pd(_mm256_mul_pd(_mm256_set1_pd(ei), ej));

                        let sigma2 = _mm256_mul_pd(sigma, sigma);
                        let sr2 = _mm256_mul_pd(sigma2, inv_r2);
                        let sr4 = _mm256_mul_pd(sr2, sr2);
                        let sr6 = _mm256_mul_pd(sr4, sr2);
                        let sr12 = _mm256_mul_pd(sr6, sr6);

                        let term = _mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(2.0), sr12), sr6);
                        let pref = _mm256_mul_pd(
                            _mm256_mul_pd(_mm256_set1_pd(24.0), epsilon),
                            _mm256_mul_pd(term, inv_r2),
                        );

                        let fx_v = _mm256_mul_pd(dx, pref);
                        let fy_v = _mm256_mul_pd(dy, pref);
                        let fz_v = _mm256_mul_pd(dz, pref);

                        let mut fx_arr = [0.0_f64; 4];
                        let mut fy_arr = [0.0_f64; 4];
                        let mut fz_arr = [0.0_f64; 4];
                        let mut r2_arr = [0.0_f64; 4];
                        _mm256_storeu_pd(fx_arr.as_mut_ptr(), fx_v);
                        _mm256_storeu_pd(fy_arr.as_mut_ptr(), fy_v);
                        _mm256_storeu_pd(fz_arr.as_mut_ptr(), fz_v);
                        _mm256_storeu_pd(r2_arr.as_mut_ptr(), r2);

                        for lane in 0..4 {
                            let jj = j + lane;
                            if jj == i || r2_arr[lane] <= 1e-24 {
                                continue;
                            }
                            fx -= fx_arr[lane];
                            fy -= fy_arr[lane];
                            fz -= fz_arr[lane];
                        }
                    }

                    j += 4;
                }

                for jj in j..n {
                    if jj == i {
                        continue;
                    }
                    let dx = xs_ref[jj] - xi;
                    let dy = ys_ref[jj] - yi;
                    let dz = zs_ref[jj] - zi;
                    let r2 = dx * dx + dy * dy + dz * dz;
                    let sigma = 0.5 * (si + sigmas_ref[jj]);
                    let epsilon = (ei * eps_ref[jj]).sqrt();
                    let pref = lj_force_prefactor_from_r2(r2, sigma, epsilon);
                    fx -= dx * pref;
                    fy -= dy * pref;
                    fz -= dz * pref;
                }

                out[local_i] = Vector3::new(fx, fy, fz);
            }

            (start, out)
        }));
    }

    let mut forces = vec![Vector3::zeros(); n];
    for handle in handles {
        let (start, block) = handle.join().expect("force worker thread panicked");
        for (k, f) in block.into_iter().enumerate() {
            forces[start + k] = f;
        }
    }

    forces
}

#[cfg(not(target_arch = "x86_64"))]
fn compute_forces_simd_parallel_avx(particles: &[Particle]) -> Vec<Vector3<f64>> {
    compute_forces_parallel_scalar(particles)
}

pub fn apply_forces_simd_parallel(particles: &mut [Particle]) {
    let forces = compute_forces_simd_parallel(particles);
    for (p, f) in particles.iter_mut().zip(forces.into_iter()) {
        p.force = f;
    }
}
