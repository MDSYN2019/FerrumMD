use crate::lennard_jones_simulations::{LJParameters, Particle};
use crate::molecule::molecule::{
    Angle, Bond, Dihedral, NonbondedPairScaling, RigidWaterDefinition, System,
};
use nalgebra::Vector3;
use std::collections::{HashMap, HashSet};
use std::fs;

#[derive(Clone, Debug, Default)]
pub struct ItpAtomType {
    pub name: String,
    pub mass: f64,
    pub charge: f64,
    pub sigma: Option<f64>,   // This is some or none
    pub epsilon: Option<f64>, // This is some or none
    pub c6: Option<f64>,      // This is some or none
    pub c12: Option<f64>,     // This is some or none
}

#[derive(Clone, Debug)]
pub struct ItpAtom {
    pub index: usize,
    pub type_name: String,
    pub charge: f64,
    pub mass: Option<f64>,
}

#[derive(Clone, Debug)]
pub struct ItpBond {
    pub atom1: usize,
    pub atom2: usize,
    pub length: f64,
    pub force_constant: f64,
}

#[derive(Clone, Debug)]
pub struct ItpAngle {
    pub atom1: usize,
    pub atom2: usize,
    pub atom3: usize,
    pub theta0_deg: f64,
    pub force_constant: f64,
}

#[derive(Clone, Debug)]
pub struct ItpDihedral {
    pub atom1: usize,
    pub atom2: usize,
    pub atom3: usize,
    pub atom4: usize,
    pub phase_deg: f64,
    pub force_constant: f64,
    pub multiplicity: usize,
}

#[derive(Clone, Debug)]
pub struct ItpSettles {
    pub oxygen_atom: usize,
    pub oh_distance: f64,
    pub hh_distance: f64,
}

#[derive(Clone, Debug)]
pub struct ItpPair {
    pub atom1: usize,
    pub atom2: usize,
    pub lj_scale: f64,
    pub coulomb_scale: f64,
}

#[derive(Clone, Debug, Default)]
pub struct ItpForceField {
    pub molecule_name: Option<String>,
    pub atom_types: HashMap<String, ItpAtomType>,
    pub atoms: Vec<ItpAtom>,
    pub bonds: Vec<ItpBond>,
    pub angles: Vec<ItpAngle>,
    pub dihedrals: Vec<ItpDihedral>,
    pub exclusions: Vec<(usize, usize)>,
    pub pairs: Vec<ItpPair>,
    pub settles: Vec<ItpSettles>,
}

impl ItpForceField {
    pub fn parse_str(contents: &str) -> Result<Self, String> {
        let mut ff = ItpForceField::default();
        let mut section = String::new();
        let defines = default_preprocessor_defines();

        for (line_number, line) in preprocess_itp_lines(contents, &defines) {
            let line = line.as_str();

            if line.starts_with('[') && line.ends_with(']') {
                section = line[1..line.len() - 1].trim().to_ascii_lowercase();
                continue;
            }

            let tokens: Vec<&str> = line.split_whitespace().collect();
            if tokens.is_empty() {
                continue;
            }

            match section.as_str() {
                "moleculetype" => {
                    if ff.molecule_name.is_none() {
                        ff.molecule_name = Some(tokens[0].to_string());
                    }
                }
                "atomtypes" => {
                    let atom_type = parse_atomtype(&tokens)
                        .map_err(|e| format!("line {}: {e}", line_number + 1))?;
                    ff.atom_types.insert(atom_type.name.clone(), atom_type);
                }
                "atoms" => {
                    ff.atoms.push(
                        parse_atom(&tokens)
                            .map_err(|e| format!("line {}: {e}", line_number + 1))?,
                    );
                }
                "bonds" => {
                    ff.bonds.push(
                        parse_bond(&tokens)
                            .map_err(|e| format!("line {}: {e}", line_number + 1))?,
                    );
                }
                "angles" => {
                    ff.angles.push(
                        parse_angle(&tokens)
                            .map_err(|e| format!("line {}: {e}", line_number + 1))?,
                    );
                }
                "dihedrals" => {
                    ff.dihedrals.push(
                        parse_dihedral(&tokens)
                            .map_err(|e| format!("line {}: {e}", line_number + 1))?,
                    );
                }
                "exclusions" => {
                    let entries = parse_exclusions(&tokens)
                        .map_err(|e| format!("line {}: {e}", line_number + 1))?;
                    ff.exclusions.extend(entries);
                }
                "pairs" => {
                    ff.pairs.push(
                        parse_pair(&tokens)
                            .map_err(|e| format!("line {}: {e}", line_number + 1))?,
                    );
                }
                "settles" => {
                    ff.settles.push(
                        parse_settles(&tokens)
                            .map_err(|e| format!("line {}: {e}", line_number + 1))?,
                    );
                }
                _ => {}
            }
        }

        if ff.atoms.is_empty() {
            return Err("itp input did not contain an [ atoms ] section".to_string());
        }

        Ok(ff)
    }

    pub fn read_itp(path: &str) -> Result<Self, String> {
        let contents = fs::read_to_string(path)
            .map_err(|e| format!("failed to read itp file at '{path}': {e}"))?;
        Self::parse_str(&contents)
    }

    pub fn read_atomtypes_from_itp(path: &str) -> Result<HashMap<String, ItpAtomType>, String> {
        let contents = fs::read_to_string(path)
            .map_err(|e| format!("failed to read itp file at '{path}': {e}"))?;
        Self::parse_atomtypes_str(&contents)
    }

    pub fn parse_atomtypes_str(contents: &str) -> Result<HashMap<String, ItpAtomType>, String> {
        let mut atom_types = HashMap::new();
        let mut section = String::new();
        let defines = default_preprocessor_defines();

        for (line_number, line) in preprocess_itp_lines(contents, &defines) {
            let line = line.as_str();

            if line.starts_with('[') && line.ends_with(']') {
                section = line[1..line.len() - 1].trim().to_ascii_lowercase();
                continue;
            }

            if section != "atomtypes" {
                continue;
            }

            let tokens: Vec<&str> = line.split_whitespace().collect();
            if tokens.is_empty() {
                continue;
            }

            let atom_type =
                parse_atomtype(&tokens).map_err(|e| format!("line {}: {e}", line_number + 1))?;
            atom_types.insert(atom_type.name.clone(), atom_type);
        }

        if atom_types.is_empty() {
            return Err("itp input did not contain an [ atomtypes ] section".to_string());
        }

        Ok(atom_types)
    }

    pub fn to_system(&self, coordinates: &[Vector3<f64>]) -> Result<System, String> {
        if coordinates.len() != self.atoms.len() {
            return Err(format!(
                "coordinate count mismatch: got {}, expected {}",
                coordinates.len(),
                self.atoms.len()
            ));
        }

        let mut particles = Vec::with_capacity(self.atoms.len());

        for (idx, atom) in self.atoms.iter().enumerate() {
            let atom_type = self
                .atom_types
                .get(&atom.type_name)
                .ok_or_else(|| format!("missing atom type '{}'", atom.type_name))?;

            let (sigma, epsilon) = infer_lj_parameters(atom_type)?;
            let mass = atom.mass.unwrap_or(atom_type.mass);

            particles.push(Particle {
                id: atom.index,
                position: coordinates[idx],
                velocity: Vector3::zeros(),
                force: Vector3::zeros(),
                lj_parameters: LJParameters {
                    epsilon,
                    sigma,
                    number_of_atoms: 1,
                },
                mass,
                energy: 0.0,
                atom_type: idx as f64,
                charge: atom.charge,
            });
        }

        let bonds = self
            .bonds
            .iter()
            .map(|b| Bond {
                atom1: b.atom1 - 1,
                atom2: b.atom2 - 1,
                k: b.force_constant,
                r0: b.length,
            })
            .collect();

        let angles = self
            .angles
            .iter()
            .map(|a| Angle {
                atom1: a.atom1 - 1,
                atom2: a.atom2 - 1,
                atom3: a.atom3 - 1,
                k: a.force_constant,
                theta0: a.theta0_deg.to_radians(),
            })
            .collect();

        let dihedrals = self
            .dihedrals
            .iter()
            .map(|d| Dihedral {
                atom1: d.atom1 - 1,
                atom2: d.atom2 - 1,
                atom3: d.atom3 - 1,
                atom4: d.atom4 - 1,
                k: d.force_constant,
                multiplicity: d.multiplicity,
                phase: d.phase_deg.to_radians(),
            })
            .collect();

        let mut system = System {
            atoms: particles,
            bonds,
            angles,
            dihedrals,
            impropers: Vec::new(),
            ..System::default()
        };

        for &(i, j) in &self.exclusions {
            system.add_exclusion(i - 1, j - 1);
        }

        for pair in &self.pairs {
            system.set_pair_scaling(
                pair.atom1 - 1,
                pair.atom2 - 1,
                NonbondedPairScaling {
                    lj_scale: pair.lj_scale,
                    coulomb_scale: pair.coulomb_scale,
                },
            );
        }

        if let Some(settle) = self.settles.first() {
            let oxygen = settle.oxygen_atom - 1;
            let hydrogens: Vec<usize> = self
                .bonds
                .iter()
                .filter_map(|bond| {
                    if bond.atom1 - 1 == oxygen {
                        Some(bond.atom2 - 1)
                    } else if bond.atom2 - 1 == oxygen {
                        Some(bond.atom1 - 1)
                    } else {
                        None
                    }
                })
                .collect();
            if hydrogens.len() == 2 {
                system.rigid_water = Some(RigidWaterDefinition {
                    oxygen_index: oxygen,
                    hydrogen1_index: hydrogens[0],
                    hydrogen2_index: hydrogens[1],
                    oh_distance: settle.oh_distance,
                    hh_distance: settle.hh_distance,
                });
            } else {
                return Err(
                    "invalid [ settles ] entry: expected oxygen bonded to two hydrogens"
                        .to_string(),
                );
            }
        }

        Ok(system)
    }
}

fn parse_atomtype(tokens: &[&str]) -> Result<ItpAtomType, String> {
    if tokens.len() < 6 {
        return Err("atomtypes row requires at least 6 columns".to_string());
    }

    let name = tokens[0].to_string();
    let ptype_index = tokens
        .iter()
        .enumerate()
        .skip(1)
        .find(|(_, token)| token.chars().all(|c| c.is_ascii_alphabetic()))
        .map(|(idx, _)| idx);

    let (mass, charge) = if let Some(ptype_idx) = ptype_index {
        if ptype_idx < 3 {
            return Err("atomtypes row missing mass/charge before ptype".to_string());
        }
        (
            parse_f64(tokens, ptype_idx - 2, "atomtype mass")?,
            parse_f64(tokens, ptype_idx - 1, "atomtype charge")?,
        )
    } else {
        let numeric_indices: Vec<usize> = tokens
            .iter()
            .enumerate()
            .skip(1)
            .filter_map(|(idx, token)| token.parse::<f64>().ok().map(|_| idx))
            .collect();
        if numeric_indices.len() < 4 {
            return Err("atomtypes row missing numeric mass/charge parameters".to_string());
        }
        (
            parse_f64(tokens, numeric_indices[0], "atomtype mass")?,
            parse_f64(tokens, numeric_indices[1], "atomtype charge")?,
        )
    };

    // Martini files can encode either sigma/epsilon or C6/C12 in the last 2 columns.
    let maybe_a = tokens
        .get(tokens.len().saturating_sub(2))
        .ok_or_else(|| "atomtype missing nonbonded columns".to_string())?
        .parse::<f64>()
        .map_err(|e| format!("invalid atomtype nonbonded value: {e}"))?;

    let maybe_b = tokens
        .get(tokens.len().saturating_sub(1))
        .ok_or_else(|| "atomtype missing nonbonded columns".to_string())?
        .parse::<f64>()
        .map_err(|e| format!("invalid atomtype nonbonded value: {e}"))?;

    // Heuristic: sigma is usually ~0.3-0.8 nm, epsilon few kJ/mol; C12 is tiny.
    let (sigma, epsilon, c6, c12) = if maybe_a > 0.0 && maybe_a < 2.0 && maybe_b < 20.0 {
        (Some(maybe_a), Some(maybe_b), None, None)
    } else {
        (None, None, Some(maybe_a), Some(maybe_b))
    };

    Ok(ItpAtomType {
        name,
        mass,
        charge,
        sigma,
        epsilon,
        c6,
        c12,
    })
}

fn parse_atom(tokens: &[&str]) -> Result<ItpAtom, String> {
    if tokens.len() < 7 {
        return Err("atoms row requires at least 7 columns".to_string());
    }

    let index = parse_usize(tokens, 0, "atom index")?;
    let type_name = tokens[1].to_string();
    let charge = parse_f64(tokens, 6, "atom charge")?;
    let mass = if tokens.len() > 7 {
        Some(parse_f64(tokens, 7, "atom mass")?)
    } else {
        None
    };

    Ok(ItpAtom {
        index,
        type_name,
        charge,
        mass,
    })
}

fn parse_bond(tokens: &[&str]) -> Result<ItpBond, String> {
    if tokens.len() < 5 {
        return Err("bonds row requires at least 5 columns".to_string());
    }

    Ok(ItpBond {
        atom1: parse_usize(tokens, 0, "bond atom1")?,
        atom2: parse_usize(tokens, 1, "bond atom2")?,
        length: parse_f64(tokens, 3, "bond length")?,
        force_constant: parse_f64(tokens, 4, "bond force constant")?,
    })
}

fn parse_angle(tokens: &[&str]) -> Result<ItpAngle, String> {
    if tokens.len() < 6 {
        return Err("angles row requires at least 6 columns".to_string());
    }

    Ok(ItpAngle {
        atom1: parse_usize(tokens, 0, "angle atom1")?,
        atom2: parse_usize(tokens, 1, "angle atom2")?,
        atom3: parse_usize(tokens, 2, "angle atom3")?,
        theta0_deg: parse_f64(tokens, 4, "angle theta0")?,
        force_constant: parse_f64(tokens, 5, "angle force constant")?,
    })
}

fn parse_dihedral(tokens: &[&str]) -> Result<ItpDihedral, String> {
    if tokens.len() < 8 {
        return Err("dihedrals row requires at least 8 columns".to_string());
    }

    Ok(ItpDihedral {
        atom1: parse_usize(tokens, 0, "dihedral atom1")?,
        atom2: parse_usize(tokens, 1, "dihedral atom2")?,
        atom3: parse_usize(tokens, 2, "dihedral atom3")?,
        atom4: parse_usize(tokens, 3, "dihedral atom4")?,
        phase_deg: parse_f64(tokens, 5, "dihedral phase")?,
        force_constant: parse_f64(tokens, 6, "dihedral force constant")?,
        multiplicity: parse_usize(tokens, 7, "dihedral multiplicity")?,
    })
}

fn parse_settles(tokens: &[&str]) -> Result<ItpSettles, String> {
    if tokens.len() < 4 {
        return Err("settles row requires at least 4 columns".to_string());
    }

    Ok(ItpSettles {
        oxygen_atom: parse_usize(tokens, 0, "settles oxygen atom")?,
        oh_distance: parse_f64(tokens, 2, "settles O-H distance")?,
        hh_distance: parse_f64(tokens, 3, "settles H-H distance")?,
    })
}

fn parse_exclusions(tokens: &[&str]) -> Result<Vec<(usize, usize)>, String> {
    if tokens.len() < 2 {
        return Err("exclusions row requires at least 2 columns".to_string());
    }
    let atom = parse_usize(tokens, 0, "exclusions atom")?;
    let mut pairs = Vec::with_capacity(tokens.len() - 1);
    for idx in 1..tokens.len() {
        let excluded_atom = parse_usize(tokens, idx, "excluded atom")?;
        pairs.push((atom, excluded_atom));
    }
    Ok(pairs)
}

fn parse_pair(tokens: &[&str]) -> Result<ItpPair, String> {
    if tokens.len() < 2 {
        return Err("pairs row requires at least 2 columns".to_string());
    }

    let lj_scale = if tokens.len() >= 5 {
        parse_f64(tokens, 3, "pair LJ scale")?
    } else {
        1.0
    };
    let coulomb_scale = if tokens.len() >= 5 {
        parse_f64(tokens, 4, "pair Coulomb scale")?
    } else {
        1.0
    };

    Ok(ItpPair {
        atom1: parse_usize(tokens, 0, "pair atom1")?,
        atom2: parse_usize(tokens, 1, "pair atom2")?,
        lj_scale,
        coulomb_scale,
    })
}

fn infer_lj_parameters(atom_type: &ItpAtomType) -> Result<(f64, f64), String> {
    // OK(T) -> An element T was found
    if let (Some(sigma), Some(epsilon)) = (atom_type.sigma, atom_type.epsilon) {
        return Ok((sigma, epsilon)); // sigma, epsilon was found and we are returning these
    }

    if let (Some(c6), Some(c12)) = (atom_type.c6, atom_type.c12) {
        // Some force fields (e.g., TIP3P hydrogens) intentionally set both C6/C12 to zero.
        if c6 == 0.0 && c12 == 0.0 {
            return Ok((0.0, 0.0));
        }
        if c6 <= 0.0 || c12 <= 0.0 {
            return Err(format!(
                "atom type '{}' has non-positive C6/C12",
                atom_type.name
            ));
        }
        let sigma = (c12 / c6).powf(1.0 / 6.0);
        let epsilon = (c6 * c6) / (4.0 * c12);
        return Ok((sigma, epsilon));
    }

    Err(format!(
        "atom type '{}' is missing LJ parameters",
        atom_type.name
    ))
}

pub type MartiniAtomType = ItpAtomType;
pub type MartiniAtom = ItpAtom;
pub type MartiniBond = ItpBond;
pub type MartiniAngle = ItpAngle;
pub type MartiniDihedral = ItpDihedral;
pub type MartiniForceField = ItpForceField;

fn parse_f64(tokens: &[&str], index: usize, label: &str) -> Result<f64, String> {
    tokens
        .get(index)
        .ok_or_else(|| format!("missing {label}"))?
        .parse::<f64>()
        .map_err(|e| format!("failed to parse {label}: {e}"))
}

fn parse_usize(tokens: &[&str], index: usize, label: &str) -> Result<usize, String> {
    tokens
        .get(index)
        .ok_or_else(|| format!("missing {label}"))?
        .parse::<usize>()
        .map_err(|e| format!("failed to parse {label}: {e}"))
}

#[derive(Clone, Copy, Debug)]
struct ConditionalFrame {
    parent_active: bool,
    condition_active: bool,
    in_else_branch: bool,
}

fn default_preprocessor_defines() -> HashSet<&'static str> {
    HashSet::from(["_FF_CHARMM", "FLEXIBLE", "CHARMM_TIP3P"])
}

fn is_active(stack: &[ConditionalFrame]) -> bool {
    stack.last().map_or(true, |frame| {
        if frame.in_else_branch {
            frame.parent_active && !frame.condition_active
        } else {
            frame.parent_active && frame.condition_active
        }
    })
}

fn preprocess_itp_lines(contents: &str, defines: &HashSet<&'static str>) -> Vec<(usize, String)> {
    let mut filtered = Vec::new();
    let mut condition_stack: Vec<ConditionalFrame> = Vec::new();

    for (line_number, raw_line) in contents.lines().enumerate() {
        let line = raw_line.split(';').next().unwrap_or("").trim();
        if line.is_empty() {
            continue;
        }

        if let Some(symbol) = line.strip_prefix("#ifdef").map(str::trim) {
            let frame = ConditionalFrame {
                parent_active: is_active(&condition_stack),
                condition_active: defines.contains(symbol),
                in_else_branch: false,
            };
            condition_stack.push(frame);
            continue;
        }

        if let Some(symbol) = line.strip_prefix("#ifndef").map(str::trim) {
            let frame = ConditionalFrame {
                parent_active: is_active(&condition_stack),
                condition_active: !defines.contains(symbol),
                in_else_branch: false,
            };
            condition_stack.push(frame);
            continue;
        }

        if line.starts_with("#else") {
            if let Some(last) = condition_stack.last_mut() {
                last.in_else_branch = true;
            }
            continue;
        }

        if line.starts_with("#endif") {
            let _ = condition_stack.pop();
            continue;
        }

        if line.starts_with('#') {
            continue;
        }

        if !is_active(&condition_stack) {
            continue;
        }

        filtered.push((line_number, line.to_string()));
    }

    filtered
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_martini_itp_and_builds_system() {
        let itp = r#"
[ moleculetype ]
; name nrexcl
WAT 1

[ atomtypes ]
;name mass charge ptype sigma epsilon
P4   72.0 0.0 A 0.47 5.0

[ atoms ]
;nr type resnr residue atom cgnr charge mass
1 P4 1 WAT W1 1 0.0 72.0
2 P4 1 WAT W2 1 0.0 72.0

[ bonds ]
1 2 1 0.30 1250

[ angles ]
1 2 1 2 120.0 25.0

[ dihedrals ]
1 2 1 2 1 180.0 5.0 1
"#;

        let ff = MartiniForceField::parse_str(itp).expect("martini parsing should succeed");
        assert_eq!(ff.atoms.len(), 2);
        assert_eq!(ff.bonds.len(), 1);

        let coords = vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(0.3, 0.0, 0.0)];
        let system = ff.to_system(&coords).expect("system build should succeed");

        assert_eq!(system.atoms.len(), 2);
        assert_eq!(system.bonds[0].atom1, 0);
        assert!((system.atoms[0].lj_parameters.sigma - 0.47).abs() < 1e-12);
        assert!((system.atoms[0].lj_parameters.epsilon - 5.0).abs() < 1e-12);
    }

    #[test]
    fn parses_c6_c12_atomtypes() {
        let itp = r#"
[ atomtypes ]
C1   72.0 0.0 A 2.0 1.0
[ atoms ]
1 C1 1 TST C1 1 0.0 72.0
"#;

        let ff = MartiniForceField::parse_str(itp).expect("martini parsing should succeed");
        let atom_type = ff.atom_types.get("C1").expect("atomtype should exist");
        let (sigma, epsilon) = infer_lj_parameters(atom_type).expect("lj inference should succeed");

        assert!(sigma > 0.0);
        assert!(epsilon > 0.0);
    }

    #[test]
    fn parses_general_gromacs_atomtypes_layout() {
        let itp = r#"
[ atomtypes ]
; name at.num mass charge ptype sigma epsilon
opls_001 6 12.011 0.0 A 0.355 0.276144
[ atoms ]
1 opls_001 1 ETH C1 1 0.0 12.011
"#;

        let ff = ItpForceField::parse_str(itp).expect("itp parsing should succeed");
        let atom_type = ff
            .atom_types
            .get("opls_001")
            .expect("atomtype should exist");

        assert!((atom_type.mass - 12.011).abs() < 1e-12);
        assert!((atom_type.charge - 0.0).abs() < 1e-12);
        let (sigma, epsilon) = infer_lj_parameters(atom_type).expect("lj inference should succeed");
        assert!((sigma - 0.355).abs() < 1e-12);
        assert!((epsilon - 0.276144).abs() < 1e-12);
    }

    #[test]
    fn parses_atomtypes_only_itp() {
        let itp = r#"
[ atomtypes ]
; name at.num mass charge ptype sigma epsilon
OWT3 8 15.999400 -0.834 A 0.315058 0.636386
"#;

        let atom_types =
            ItpForceField::parse_atomtypes_str(itp).expect("atomtypes parsing should succeed");
        let atom_type = atom_types.get("OWT3").expect("OWT3 should exist");
        assert!((atom_type.mass - 15.9994).abs() < 1e-12);
        assert!((atom_type.charge + 0.834).abs() < 1e-12);
    }

    #[test]
    fn parses_settles_exclusions_and_pairs() {
        let itp = r#"
[ atomtypes ]
OW 15.9994 -0.834 A 0.315058 0.636386
HW 1.008 0.417 A 0.0 0.0
[ atoms ]
1 OW 1 WAT OW 1 -0.834 15.9994
2 HW 1 WAT HW1 1 0.417 1.008
3 HW 1 WAT HW2 1 0.417 1.008
[ bonds ]
1 2 1 0.09572 1000
1 3 1 0.09572 1000
[ settles ]
1 1 0.09572 0.15139
[ exclusions ]
1 2 3
[ pairs ]
2 3 1 0.5 0.8333
"#;
        let ff = ItpForceField::parse_str(itp).expect("itp parsing should succeed");
        assert_eq!(ff.settles.len(), 1);
        assert_eq!(ff.exclusions.len(), 2);
        assert_eq!(ff.pairs.len(), 1);
        let coords = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(0.09572, 0.0, 0.0),
            Vector3::new(-0.023999, 0.092663, 0.0),
        ];
        let system = ff.to_system(&coords).expect("system build should succeed");
        assert!(system.rigid_water.is_some());
        assert!(system.exclusions.contains(&(0, 1)));
        assert!(system.exclusions.contains(&(0, 2)));
        let pair = system
            .pair_scalings
            .get(&(1, 2))
            .expect("pair scaling should exist");
        assert!((pair.lj_scale - 0.5).abs() < 1e-12);
        assert!((pair.coulomb_scale - 0.8333).abs() < 1e-12);
    }
}
