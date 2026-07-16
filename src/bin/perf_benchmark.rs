use nalgebra::Vector3;
use rand::Rng;
use sang_md::lennard_jones_simulations::{LJParameters, Particle};
use sang_md::performance::{compute_forces_python_baseline, compute_forces_simd_parallel};
use std::time::Instant;

fn make_particles(n: usize, box_length: f64) -> Vec<Particle> {
    let mut rng = rand::rng();
    (0..n)
        .map(|id| Particle {
            id,
            position: Vector3::new(
                rng.random_range(0.0..box_length),
                rng.random_range(0.0..box_length),
                rng.random_range(0.0..box_length),
            ),
            velocity: Vector3::zeros(),
            force: Vector3::zeros(),
            lj_parameters: LJParameters {
                epsilon: 1.0,
                sigma: 1.0,
                number_of_atoms: 1,
            },
            mass: 1.0,
            energy: 0.0,
            atom_type: 1.0,
            charge: 0.0,
        })
        .collect()
}

fn avg_ms<F: Fn()>(iters: usize, f: F) -> f64 {
    let t0 = Instant::now();
    for _ in 0..iters {
        f();
    }
    t0.elapsed().as_secs_f64() * 1_000.0 / (iters as f64)
}

fn main() {
    let n = 1024usize;
    let iters = 5usize;
    let particles = make_particles(n, 20.0);

    let baseline_ms = avg_ms(iters, || {
        let _ = compute_forces_python_baseline(&particles);
    });

    let optimized_ms = avg_ms(iters, || {
        let _ = compute_forces_simd_parallel(&particles);
    });

    let speedup = baseline_ms / optimized_ms;

    println!("Particle count: {n}");
    println!(
        "Baseline (Python-style nested-loop) : {:.3} ms",
        baseline_ms
    );
    println!(
        "Optimized (SIMD + multi-threading)            : {:.3} ms",
        optimized_ms
    );
    println!("Speedup over baseline               : {:.2}x", speedup);
    println!(
        "OpenMM comparison: supply your local OpenMM timing and divide by optimized timing to compute relative speedup."
    );
}
