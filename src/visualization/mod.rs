use crate::lennard_jones_simulations::Particle;
use crate::molecule::io::read_gro;
use nalgebra::Vector3;
use xdrfile::{Frame, Trajectory, XTCTrajectory};

#[derive(Clone, Debug)]
pub struct TrajectoryFrame {
    pub step: usize,
    pub time_ps: f32,
    pub positions_nm: Vec<Vector3<f32>>,
}

#[derive(Clone, Debug)]
pub struct TrajectoryData {
    pub atom_count: usize,
    pub box_dims_nm: Option<Vector3<f32>>,
    pub frames: Vec<TrajectoryFrame>,
}

fn particles_to_positions(particles: &[Particle]) -> Vec<Vector3<f32>> {
    particles
        .iter()
        .map(|particle| {
            Vector3::new(
                particle.position.x as f32,
                particle.position.y as f32,
                particle.position.z as f32,
            )
        })
        .collect()
}

pub fn load_trajectory(gro_path: &str, xtc_path: Option<&str>) -> Result<TrajectoryData, String> {
    let (particles, gro_box_dims) = read_gro(gro_path)?;
    let atom_count = particles.len();

    if atom_count == 0 {
        return Err("GRO file has zero atoms; cannot visualize an empty system".to_string());
    }

    let mut frames = Vec::new();
    let mut box_dims_nm = gro_box_dims.map(|v| Vector3::new(v.x as f32, v.y as f32, v.z as f32));

    if let Some(path) = xtc_path {
        let mut trajectory =
            XTCTrajectory::open_read(path).map_err(|e| format!("failed to open xtc file: {e}"))?;

        let xtc_atom_count = trajectory
            .get_num_atoms()
            .map_err(|e| format!("failed to read atom count from xtc: {e}"))?;

        if xtc_atom_count != atom_count {
            return Err(format!(
                "atom count mismatch between GRO ({atom_count}) and XTC ({xtc_atom_count})"
            ));
        }

        loop {
            let mut frame = Frame::with_len(atom_count);
            match trajectory.read(&mut frame) {
                Ok(()) => {
                    if box_dims_nm.is_none() {
                        box_dims_nm = Some(Vector3::new(
                            frame.box_vector[0][0],
                            frame.box_vector[1][1],
                            frame.box_vector[2][2],
                        ));
                    }

                    frames.push(TrajectoryFrame {
                        step: frame.step,
                        time_ps: frame.time,
                        positions_nm: frame
                            .coords
                            .iter()
                            .map(|coord| Vector3::new(coord[0], coord[1], coord[2]))
                            .collect(),
                    });
                }
                Err(err) if err.is_eof() => break,
                Err(err) => {
                    return Err(format!("failed while reading xtc trajectory frame: {err}"));
                }
            }
        }
    }

    if frames.is_empty() {
        frames.push(TrajectoryFrame {
            step: 0,
            time_ps: 0.0,
            positions_nm: particles_to_positions(&particles),
        });
    }

    Ok(TrajectoryData {
        atom_count,
        box_dims_nm,
        frames,
    })
}
