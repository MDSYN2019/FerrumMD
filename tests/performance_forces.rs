use nalgebra::Vector3;
use sang_md::lennard_jones_simulations::{LJParameters, Particle};
use sang_md::performance::{compute_forces_python_baseline, compute_forces_simd_parallel};

fn particle(id: usize, x: f64, y: f64, z: f64) -> Particle {
    Particle {
        id,
        position: Vector3::new(x, y, z),
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
    }
}

#[test]
fn simd_parallel_matches_baseline_forces() {
    let particles = vec![
        particle(0, 0.2, 0.3, 0.4),
        particle(1, 1.1, 0.8, 0.6),
        particle(2, 0.7, 1.4, 1.2),
        particle(3, 1.8, 1.1, 0.9),
        particle(4, 2.2, 0.4, 1.7),
    ];

    let baseline = compute_forces_python_baseline(&particles);
    let optimized = compute_forces_simd_parallel(&particles);

    for (fb, fo) in baseline.iter().zip(optimized.iter()) {
        let diff = (fb - fo).norm();
        assert!(
            diff < 1e-8,
            "force mismatch: baseline={fb:?}, optimized={fo:?}, diff={diff}"
        );
    }
}
