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

        // GRO coordinates are in nm and FerrumMD stores coordinates in nm.
        let position = Vector3::new(x_nm, y_nm, z_nm);

        particles.push(particle_from_coordinates(atom_idx + 1, atom_name, position));
    }

    let box_line = lines[2 + natoms].trim();
    let box_values: Vec<f64> = box_line
        .split_whitespace()
        .map(|v| v.parse::<f64>())
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| format!("failed to parse gro box line: {e}"))?;

    let box_dims = if box_values.len() >= 3 {
        Some(Vector3::new(box_values[0], box_values[1], box_values[2]))
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

fn parse_lammps_axis_bounds(line: &str, axis: &str) -> Result<(f64, f64), String> {
    let fields: Vec<&str> = line.split_whitespace().collect();
    if fields.len() < 2 {
        return Err(format!("invalid {axis} bounds line in lammps data: {line}"));
    }

    let lo = fields[0]
        .parse::<f64>()
        .map_err(|e| format!("failed to parse {axis}lo bound in lammps data: {e}; line: {line}"))?;
    let hi = fields[1]
        .parse::<f64>()
        .map_err(|e| format!("failed to parse {axis}hi bound in lammps data: {e}; line: {line}"))?;
    Ok((lo, hi))
}

fn parse_lammps_atoms_line(line: &str, id: usize) -> Result<Particle, String> {
    let fields: Vec<&str> = line.split_whitespace().collect();
    if fields.len() < 4 {
        return Err(format!("invalid lammps atom line (expected x y z fields): {line}"));
    }

    // Typical "Atoms" formats are either:
    // 1) atom-ID atom-type x y z
    // 2) atom-ID mol-ID atom-type q x y z
    let coordinate_start = if fields.len() >= 7 { 4 } else { fields.len() - 3 };
    let x = fields[coordinate_start]
        .parse::<f64>()
        .map_err(|e| format!("failed to parse x coordinate in lammps data atom line: {e}; line: {line}"))?;
    let y = fields[coordinate_start + 1]
        .parse::<f64>()
        .map_err(|e| format!("failed to parse y coordinate in lammps data atom line: {e}; line: {line}"))?;
    let z = fields[coordinate_start + 2]
        .parse::<f64>()
        .map_err(|e| format!("failed to parse z coordinate in lammps data atom line: {e}; line: {line}"))?;

    Ok(particle_from_coordinates(id, "X", Vector3::new(x, y, z)))
}

pub fn read_lammps_data_from_str(
    contents: &str,
) -> Result<(Vec<Particle>, Option<Vector3<f64>>), String> {
    let lines: Vec<&str> = contents.lines().collect();

    let mut x_bounds: Option<(f64, f64)> = None;
    let mut y_bounds: Option<(f64, f64)> = None;
    let mut z_bounds: Option<(f64, f64)> = None;

    for line in &lines {
        if line.contains("xlo") && line.contains("xhi") {
            x_bounds = Some(parse_lammps_axis_bounds(line, "x")?);
        } else if line.contains("ylo") && line.contains("yhi") {
            y_bounds = Some(parse_lammps_axis_bounds(line, "y")?);
        } else if line.contains("zlo") && line.contains("zhi") {
            z_bounds = Some(parse_lammps_axis_bounds(line, "z")?);
        }
    }

    let atoms_header_idx = lines
        .iter()
        .position(|line| line.trim().starts_with("Atoms"))
        .ok_or_else(|| "failed to find 'Atoms' section in lammps data file".to_string())?;

    let mut particles = Vec::new();
    for line in lines.iter().skip(atoms_header_idx + 1) {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        if trimmed.chars().next().is_some_and(|c| c.is_ascii_alphabetic()) {
            break;
        }

        let id = particles.len() + 1;
        particles.push(parse_lammps_atoms_line(trimmed, id)?);
    }

    if particles.is_empty() {
        return Err("lammps data file has no atoms in 'Atoms' section".to_string());
    }

    let box_dims = match (x_bounds, y_bounds, z_bounds) {
        (Some((xlo, xhi)), Some((ylo, yhi)), Some((zlo, zhi))) => {
            Some(Vector3::new(xhi - xlo, yhi - ylo, zhi - zlo))
        }
        _ => None,
    };

    Ok((particles, box_dims))
}

pub fn read_lammps_data(path: &str) -> Result<(Vec<Particle>, Option<Vector3<f64>>), String> {
    let contents = fs::read_to_string(path)
        .map_err(|e| format!("failed to read lammps data file at '{path}': {e}"))?;
    read_lammps_data_from_str(&contents)
}

pub fn read_lammps_dump_from_str(
    contents: &str,
) -> Result<(Vec<Particle>, Option<Vector3<f64>>), String> {
    let lines: Vec<&str> = contents.lines().collect();
    let bounds_header_idx = lines
        .iter()
        .position(|line| line.trim() == "ITEM: BOX BOUNDS pp pp pp")
        .or_else(|| {
            lines
                .iter()
                .position(|line| line.trim().starts_with("ITEM: BOX BOUNDS"))
        });

    let mut box_dims = None;
    if let Some(idx) = bounds_header_idx {
        if lines.len() >= idx + 4 {
            let (xlo, xhi) = parse_lammps_axis_bounds(lines[idx + 1], "x")?;
            let (ylo, yhi) = parse_lammps_axis_bounds(lines[idx + 2], "y")?;
            let (zlo, zhi) = parse_lammps_axis_bounds(lines[idx + 3], "z")?;
            box_dims = Some(Vector3::new(xhi - xlo, yhi - ylo, zhi - zlo));
        }
    }

    let atoms_header_idx = lines
        .iter()
        .position(|line| line.trim().starts_with("ITEM: ATOMS"))
        .ok_or_else(|| "failed to find 'ITEM: ATOMS' section in lammps dump".to_string())?;

    let header_fields: Vec<&str> = lines[atoms_header_idx]
        .split_whitespace()
        .skip(2)
        .collect();
    let x_idx = header_fields
        .iter()
        .position(|f| *f == "x" || *f == "xu" || *f == "xs")
        .ok_or_else(|| "lammps dump ATOMS header missing x/xu/xs column".to_string())?;
    let y_idx = header_fields
        .iter()
        .position(|f| *f == "y" || *f == "yu" || *f == "ys")
        .ok_or_else(|| "lammps dump ATOMS header missing y/yu/ys column".to_string())?;
    let z_idx = header_fields
        .iter()
        .position(|f| *f == "z" || *f == "zu" || *f == "zs")
        .ok_or_else(|| "lammps dump ATOMS header missing z/zu/zs column".to_string())?;

    let mut particles = Vec::new();
    for line in lines.iter().skip(atoms_header_idx + 1) {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if trimmed.starts_with("ITEM:") {
            break;
        }

        let fields: Vec<&str> = trimmed.split_whitespace().collect();
        let max_index = *[x_idx, y_idx, z_idx].iter().max().expect("non-empty");
        if fields.len() <= max_index {
            return Err(format!("invalid lammps dump atom line: {trimmed}"));
        }

        let x = fields[x_idx]
            .parse::<f64>()
            .map_err(|e| format!("failed to parse x in lammps dump atom line: {e}; line: {trimmed}"))?;
        let y = fields[y_idx]
            .parse::<f64>()
            .map_err(|e| format!("failed to parse y in lammps dump atom line: {e}; line: {trimmed}"))?;
        let z = fields[z_idx]
            .parse::<f64>()
            .map_err(|e| format!("failed to parse z in lammps dump atom line: {e}; line: {trimmed}"))?;

        particles.push(particle_from_coordinates(
            particles.len() + 1,
            "X",
            Vector3::new(x, y, z),
        ));
    }

    if particles.is_empty() {
        return Err("lammps dump has no atoms in ITEM: ATOMS section".to_string());
    }

    Ok((particles, box_dims))
}

pub fn read_lammps_dump(path: &str) -> Result<(Vec<Particle>, Option<Vector3<f64>>), String> {
    let contents = fs::read_to_string(path)
        .map_err(|e| format!("failed to read lammps dump file at '{path}': {e}"))?;
    read_lammps_dump_from_str(&contents)
}

pub fn write_gro(
    path: &str,
    particles: &[Particle],
    box_dims: Vector3<f64>,
    title: &str,
) -> Result<(), String> {
    // GRO expects distances in nm and FerrumMD stores positions in nm.
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
            particle.position.x,
            particle.position.y,
            particle.position.z,
        ));
    }

    output.push_str(&format!(
        "{:>10.5}{:>10.5}{:>10.5}\n",
        box_dims.x, box_dims.y, box_dims.z,
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
                atom.position.x,
                atom.position.y,
                atom.position.z,
            ));

            atom_serial += 1;
        }
    }

    output.push_str(&format!(
        "{:>10.5}{:>10.5}{:>10.5}\n",
        box_dims.x, box_dims.y, box_dims.z,
    ));

    fs::write(path, output).map_err(|e| format!("failed to write gro file at '{path}': {e}"))
}

pub fn systems_to_particles_frame(systems: &[System]) -> Vec<Particle> {
    let natoms: usize = systems.iter().map(|system| system.atoms.len()).sum();
    let mut frame = Vec::with_capacity(natoms);
    for system in systems {
        frame.extend(system.atoms.iter().cloned());
    }
    frame
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
        [box_dims.x as f32, 0.0, 0.0],
        [0.0, box_dims.y as f32, 0.0],
        [0.0, 0.0, box_dims.z as f32],
    ];

    for (step, frame_particles) in frames.iter().enumerate() {
        let mut frame = Frame::with_len(natoms);
        frame.step = step;
        frame.time = step as f32 * dt_ps;
        frame.box_vector = box_nm;

        for (atom_idx, particle) in frame_particles.iter().enumerate() {
            frame.coords[atom_idx] = [
                particle.position.x as f32,
                particle.position.y as f32,
                particle.position.z as f32,
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
        assert!((particles[0].position.x - 0.111).abs() < 1e-9);

        let box_dims = box_dims.expect("box dims should exist");
        assert!((box_dims.x - 1.0).abs() < 1e-9);
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
        assert!(content.contains("1.000"));
    }

    #[test]
    fn parses_lammps_data_atoms_and_box() {
        let data = "LAMMPS data file\n\
\n\
2 atoms\n\
\n\
0.0 10.0 xlo xhi\n\
0.0 20.0 ylo yhi\n\
0.0 30.0 zlo zhi\n\
\n\
Atoms # atomic\n\
\n\
1 1 1.0 2.0 3.0\n\
2 1 4.0 5.0 6.0\n\
";

        let (particles, box_dims) =
            read_lammps_data_from_str(data).expect("lammps data should parse");
        assert_eq!(particles.len(), 2);
        assert!((particles[0].position.x - 1.0).abs() < 1e-9);
        assert!((particles[1].position.z - 6.0).abs() < 1e-9);

        let box_dims = box_dims.expect("box dims should exist");
        assert!((box_dims.x - 10.0).abs() < 1e-9);
        assert!((box_dims.y - 20.0).abs() < 1e-9);
        assert!((box_dims.z - 30.0).abs() < 1e-9);
    }

    #[test]
    fn parses_lammps_dump_atoms_and_box() {
        let dump = "ITEM: TIMESTEP\n\
0\n\
ITEM: NUMBER OF ATOMS\n\
2\n\
ITEM: BOX BOUNDS pp pp pp\n\
0.0 8.0\n\
0.0 9.0\n\
0.0 10.0\n\
ITEM: ATOMS id type x y z\n\
1 1 0.1 0.2 0.3\n\
2 1 0.4 0.5 0.6\n\
";

        let (particles, box_dims) =
            read_lammps_dump_from_str(dump).expect("lammps dump should parse");
        assert_eq!(particles.len(), 2);
        assert!((particles[0].position.y - 0.2).abs() < 1e-9);
        assert!((particles[1].position.x - 0.4).abs() < 1e-9);

        let box_dims = box_dims.expect("box dims should exist");
        assert!((box_dims.x - 8.0).abs() < 1e-9);
        assert!((box_dims.y - 9.0).abs() < 1e-9);
        assert!((box_dims.z - 10.0).abs() < 1e-9);
    }
}
