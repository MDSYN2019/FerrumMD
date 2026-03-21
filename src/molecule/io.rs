use crate::lennard_jones_simulations::{LJParameters, Particle};
use crate::molecule::molecule::System;
use nalgebra::Vector3;
use std::fs;
use xdrfile::{Frame, Trajectory, XTCTrajectory};

fn default_mass(atom_name: &str) -> f64 {
    let trimmed = atom_name.trim();
    let symbol = trimmed
        .chars()
        .find(|c| c.is_ascii_alphabetic())
        .map(|c| c.to_ascii_uppercase());

    match symbol {
        Some('H') => 1.008,
        Some('C') => 12.011,
        Some('N') => 14.007,
        Some('O') => 15.999,
        Some('P') => 30.974,
        Some('S') => 32.06,
        _ => 12.0,
    }
}

fn particle_from_coordinates(id: usize, atom_name: &str, position: Vector3<f64>) -> Particle {
    Particle {
        id,
        position,
        velocity: Vector3::zeros(),
        force: Vector3::zeros(),
        lj_parameters: LJParameters {
            epsilon: 1.0,
            sigma: 1.0,
            number_of_atoms: 1,
        },
        mass: default_mass(atom_name),
        energy: 0.0,
        atom_type: 0.0,
        charge: 0.0,
    }
}

fn parse_f64_slice(line: &str, start: usize, end: usize, label: &str) -> Result<f64, String> {
    line.get(start..end)
        .ok_or_else(|| format!("missing {label} field in line: {line}"))?
        .trim()
        .parse::<f64>()
        .map_err(|e| format!("failed to parse {label}: {e}; line: {line}"))
}

pub fn read_pdb_from_str(contents: &str) -> Result<Vec<Particle>, String> {
    let mut particles = Vec::new();

    for line in contents.lines() {
        if line.starts_with("ATOM") || line.starts_with("HETATM") {
            let atom_name = line.get(12..16).unwrap_or("X").trim();
            let x = parse_f64_slice(line, 30, 38, "x")?;
            let y = parse_f64_slice(line, 38, 46, "y")?;
            let z = parse_f64_slice(line, 46, 54, "z")?;

            particles.push(particle_from_coordinates(
                particles.len() + 1,
                atom_name,
                Vector3::new(x, y, z),
            ));
        }
    }

    if particles.is_empty() {
        return Err("no ATOM/HETATM records found in pdb input".to_string());
    }

    Ok(particles)
}

pub fn read_pdb(path: &str) -> Result<Vec<Particle>, String> {
    let contents = fs::read_to_string(path)
        .map_err(|e| format!("failed to read pdb file at '{path}': {e}"))?;
    read_pdb_from_str(&contents)
}

pub fn read_gro_from_str(contents: &str) -> Result<(Vec<Particle>, Option<Vector3<f64>>), String> {
    let lines: Vec<&str> = contents.lines().collect();
    if lines.len() < 3 {
        return Err("gro input must contain at least 3 lines".to_string());
    }

    let natoms = lines[1]
        .trim()
        .parse::<usize>()
        .map_err(|e| format!("failed to parse atom count from gro header: {e}"))?;

    if lines.len() < natoms + 3 {
        return Err(format!(
            "gro input declares {natoms} atoms but only {} lines provided",
            lines.len()
        ));
    }

    let mut particles = Vec::with_capacity(natoms);
    for atom_idx in 0..natoms {
        let line = lines[2 + atom_idx];

        let atom_name = line.get(10..15).unwrap_or("X").trim();
        let parsed_fixed_x = parse_f64_slice(line, 20, 28, "x");
        let parsed_fixed_y = parse_f64_slice(line, 28, 36, "y");
        let parsed_fixed_z = parse_f64_slice(line, 36, 44, "z");

        let (x_nm, y_nm, z_nm) = match (parsed_fixed_x, parsed_fixed_y, parsed_fixed_z) {
            (Ok(x), Ok(y), Ok(z)) => (x, y, z),
            _ => {
                let tokens: Vec<&str> = line.split_whitespace().collect();
                if tokens.len() < 3 {
                    return Err(format!(
                        "failed to parse coordinates from gro atom line: {line}"
                    ));
                }
                let n = tokens.len();
                let x = tokens[n - 3]
                    .parse::<f64>()
                    .map_err(|e| format!("failed to parse x: {e}; line: {line}"))?;
                let y = tokens[n - 2]
                    .parse::<f64>()
                    .map_err(|e| format!("failed to parse y: {e}; line: {line}"))?;
                let z = tokens[n - 1]
                    .parse::<f64>()
                    .map_err(|e| format!("failed to parse z: {e}; line: {line}"))?;
                (x, y, z)
            }
        };

        // GRO coordinates are in nm and the rest of this repository uses reduced/Å-like units.
        // Keep units explicit by converting nm -> Å.
        let position = Vector3::new(x_nm * 10.0, y_nm * 10.0, z_nm * 10.0);

        particles.push(particle_from_coordinates(atom_idx + 1, atom_name, position));
    }

    let box_line = lines[2 + natoms].trim();
    let box_values: Vec<f64> = box_line
        .split_whitespace()
        .map(|v| v.parse::<f64>())
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| format!("failed to parse gro box line: {e}"))?;

    let box_dims = if box_values.len() >= 3 {
        Some(Vector3::new(
            box_values[0] * 10.0,
            box_values[1] * 10.0,
            box_values[2] * 10.0,
        ))
    } else {
        None
    };

    Ok((particles, box_dims))
}

pub fn read_gro(path: &str) -> Result<(Vec<Particle>, Option<Vector3<f64>>), String> {
    let contents = fs::read_to_string(path)
        .map_err(|e| format!("failed to read gro file at '{path}': {e}"))?;
    read_gro_from_str(&contents)
}

pub fn write_gro(
    path: &str,
    particles: &[Particle],
    box_dims: Vector3<f64>,
    title: &str,
) -> Result<(), String> {
    // GRO expects distances in nm while this codebase stores positions in Å-like units.
    let mut output = String::new();
    output.push_str(title);
    output.push('\n');
    output.push_str(&format!("{:>5}\n", particles.len()));

    for (index, particle) in particles.iter().enumerate() {
        output.push_str(&format!(
            "{:>5}{:<5}{:>5}{:>5}{:>8.3}{:>8.3}{:>8.3}\n",
            1,
            "WAT",
            "OW",
            index + 1,
            particle.position.x / 10.0,
            particle.position.y / 10.0,
            particle.position.z / 10.0,
        ));
    }

    output.push_str(&format!(
        "{:>10.5}{:>10.5}{:>10.5}\n",
        box_dims.x / 10.0,
        box_dims.y / 10.0,
        box_dims.z / 10.0,
    ));

    fs::write(path, output).map_err(|e| format!("failed to write gro file at '{path}': {e}"))
}

pub fn write_gro_systems(
    path: &str,
    systems: &[System],
    box_dims: Vector3<f64>,
    title: &str,
) -> Result<(), String> {
    let natoms: usize = systems.iter().map(|system| system.atoms.len()).sum();

    let mut output = String::new();
    output.push_str(title);
    output.push('\n');
    output.push_str(&format!("{:>5}\n", natoms));

    let mut atom_serial = 1usize;
    for (res_index, system) in systems.iter().enumerate() {
        for (atom_idx, atom) in system.atoms.iter().enumerate() {
            let atom_name = match atom_idx {
                0 => "OW",
                1 => "HW1",
                2 => "HW2",
                _ => "X",
            };

            output.push_str(&format!(
                "{:>5}{:<5}{:>5}{:>5}{:>8.3}{:>8.3}{:>8.3}\n",
                (res_index + 1) % 100_000,
                "WAT",
                atom_name,
                atom_serial % 100_000,
                atom.position.x / 10.0,
                atom.position.y / 10.0,
                atom.position.z / 10.0,
            ));

            atom_serial += 1;
        }
    }

    output.push_str(&format!(
        "{:>10.5}{:>10.5}{:>10.5}\n",
        box_dims.x / 10.0,
        box_dims.y / 10.0,
        box_dims.z / 10.0,
    ));

    fs::write(path, output).map_err(|e| format!("failed to write gro file at '{path}': {e}"))
}

pub fn write_xtc(
    path: &str,
    frames: &[Vec<Particle>],
    box_dims: Vector3<f64>,
    dt_ps: f32,
) -> Result<(), String> {
    if frames.is_empty() {
        return Err("cannot write xtc with zero frames".to_string());
    }

    let natoms = frames[0].len();
    if natoms == 0 {
        return Err("cannot write xtc with zero atoms".to_string());
    }

    for (idx, frame) in frames.iter().enumerate() {
        if frame.len() != natoms {
            return Err(format!(
                "xtc frame {idx} has {} atoms but expected {natoms}",
                frame.len()
            ));
        }
    }

    let mut trajectory =
        XTCTrajectory::open_write(path).map_err(|e| format!("failed to open xtc file: {e}"))?;

    let box_nm = [
        [box_dims.x as f32 / 10.0, 0.0, 0.0],
        [0.0, box_dims.y as f32 / 10.0, 0.0],
        [0.0, 0.0, box_dims.z as f32 / 10.0],
    ];

    for (step, frame_particles) in frames.iter().enumerate() {
        let mut frame = Frame::with_len(natoms);
        frame.step = step;
        frame.time = step as f32 * dt_ps;
        frame.box_vector = box_nm;

        for (atom_idx, particle) in frame_particles.iter().enumerate() {
            frame.coords[atom_idx] = [
                particle.position.x as f32 / 10.0,
                particle.position.y as f32 / 10.0,
                particle.position.z as f32 / 10.0,
            ];
        }

        trajectory
            .write(&frame)
            .map_err(|e| format!("failed to write xtc frame {step}: {e}"))?;
    }

    trajectory
        .flush()
        .map_err(|e| format!("failed to flush xtc trajectory: {e}"))
}

pub fn write_xtc_system(
    path: &str,
    frames: &[Vec<System>],
    box_dims: Vector3<f64>,
    dt_ps: f32,
) -> Result<(), String> {
    if frames.is_empty() {
        return Err("cannot write xtc with zero frames".to_string());
    }

    let natoms = frames[0].len();
    if natoms == 0 {
        return Err("cannot write xtc with zero atoms".to_string());
    }

    for (idx, frame) in frames.iter().enumerate() {
        if frame.len() != natoms {
            return Err(format!(
                "xtc frame {idx} has {} atoms but expected {natoms}",
                frame.len()
            ));
        }
    }

    let mut trajectory =
        XTCTrajectory::open_write(path).map_err(|e| format!("failed to open xtc file: {e}"))?;

    let box_nm = [
        [box_dims.x as f32 / 10.0, 0.0, 0.0],
        [0.0, box_dims.y as f32 / 10.0, 0.0],
        [0.0, 0.0, box_dims.z as f32 / 10.0],
    ];

    for (step, frame_particles) in frames.iter().enumerate() {
        let mut frame = Frame::with_len(natoms);
        frame.step = step;
        frame.time = step as f32 * dt_ps;
        frame.box_vector = box_nm;

        for (atom_idx, particle) in frame_particles.iter().enumerate() {
            frame.coords[atom_idx] = [
                particle.position.x as f32 / 10.0,
                particle.position.y as f32 / 10.0,
                particle.position.z as f32 / 10.0,
            ];
        }

        trajectory
            .write(&frame)
            .map_err(|e| format!("failed to write xtc frame {step}: {e}"))?;
    }

    trajectory
        .flush()
        .map_err(|e| format!("failed to flush xtc trajectory: {e}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_basic_pdb_records() {
        let pdb =
            "ATOM      1  O   HOH A   1      11.104  13.207  -8.188  1.00 10.00           O\n\
ATOM      2  H1  HOH A   1      10.321  13.726  -8.003  1.00 10.00           H\n\
END\n";

        let particles = read_pdb_from_str(pdb).expect("pdb should parse");
        assert_eq!(particles.len(), 2);
        assert!((particles[0].position.x - 11.104).abs() < 1e-9);
        assert_eq!(particles[1].id, 2);
        assert!((particles[1].mass - 1.008).abs() < 1e-6);
    }

    #[test]
    fn parses_basic_gro_records() {
        let gro = "Test system\n\
2\n\
    1WAT     O    1   0.111   0.222   0.333\n\
    1WAT    H1    2   0.121   0.232   0.343\n\
   1.00000   1.00000   1.00000\n";

        let (particles, box_dims) = read_gro_from_str(gro).expect("gro should parse");
        assert_eq!(particles.len(), 2);
        assert!((particles[0].position.x - 1.11).abs() < 1e-9);

        let box_dims = box_dims.expect("box dims should exist");
        assert!((box_dims.x - 10.0).abs() < 1e-9);
    }

    #[test]
    fn writes_gro_records() {
        let particles = vec![particle_from_coordinates(
            1,
            "O",
            Vector3::new(1.0, 2.0, 3.0),
        )];

        let path = std::env::temp_dir().join(format!(
            "sang_md_test_{}_{}.gro",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .expect("time")
                .as_nanos()
        ));

        write_gro(
            path.to_str().expect("utf8 temp path"),
            &particles,
            Vector3::new(10.0, 10.0, 10.0),
            "Test",
        )
        .expect("gro write should succeed");

        let content = fs::read_to_string(&path).expect("gro file readable");
        let _ = fs::remove_file(path);
        assert!(content.contains("WAT"));
        assert!(content.contains("0.100"));
    }
}
