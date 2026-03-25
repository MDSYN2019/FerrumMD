use crate::lennard_jones_simulations::{LJParameters, Particle};
use crate::molecule::molecule::{Angle, Bond, Dihedral, Improper, System};
use nalgebra::Vector3;
use std::collections::HashMap;
use std::fs;

#[derive(Clone, Debug, Default)]
pub struct CharmmAtomType {
    pub name: String,
    pub mass: f64,
    pub element: Option<String>,
    pub epsilon: Option<f64>,
    pub sigma: Option<f64>,
    pub epsilon_14: Option<f64>,
    pub sigma_14: Option<f64>,
}

#[derive(Clone, Debug)]
pub struct CharmmAtom {
    pub index: usize,
    pub name: String,
    pub type_name: String,
    pub charge: f64,
    pub mass: Option<f64>,
}

#[derive(Clone, Debug)]
pub struct CharmmBondDef {
    pub atom1: usize,
    pub atom2: usize,
}

#[derive(Clone, Debug)]
pub struct CharmmAngleDef {
    pub atom1: usize,
    pub atom2: usize,
    pub atom3: usize,
}

#[derive(Clone, Debug)]
pub struct CharmmDihedralDef {
    pub atom1: usize,
    pub atom2: usize,
    pub atom3: usize,
    pub atom4: usize,
}

#[derive(Clone, Debug)]
struct BondParam {
    t1: String,
    t2: String,
    k: f64,
    r0: f64,
}

#[derive(Clone, Debug)]
struct AngleParam {
    t1: String,
    t2: String,
    t3: String,
    k: f64,
    theta0_deg: f64,
}

#[derive(Clone, Debug)]
struct DihedralParam {
    t1: String,
    t2: String,
    t3: String,
    t4: String,
    k: f64,
    multiplicity: usize,
    phase_deg: f64,
}

#[derive(Clone, Debug)]
struct ImproperParam {
    t1: String,
    t2: String,
    t3: String,
    t4: String,
    k: f64,
    psi0_deg: f64,
}

#[derive(Clone, Debug, Default)]
pub struct CharmmForceField {
    pub residue_name: Option<String>,
    pub atom_types: HashMap<String, CharmmAtomType>,
    pub atoms: Vec<CharmmAtom>,
    pub bonds: Vec<CharmmBondDef>,
    pub angles: Vec<CharmmAngleDef>,
    pub dihedrals: Vec<CharmmDihedralDef>,
    pub impropers: Vec<CharmmDihedralDef>,
    bond_params: Vec<BondParam>,
    angle_params: Vec<AngleParam>,
    dihedral_params: Vec<DihedralParam>,
    improper_params: Vec<ImproperParam>,
}

#[derive(Copy, Clone, Eq, PartialEq)]
enum ParamSection {
    None,
    Bonds,
    Angles,
    Dihedrals,
    Impropers,
    Nonbonded,
}

impl CharmmForceField {
    //pub fn parse_file(path_to_charmm_dir: &str) -> String {
    //    let charmm_contents = fs::read_to_string().expect("Should be able to read the charmm file");
    //}

    pub fn parse_str(contents: &str) -> Result<Self, String> {
        let mut ff = CharmmForceField::default();
        let mut atom_index_by_name: HashMap<String, usize> = HashMap::new();
        let mut section = ParamSection::None;

        for (line_number, raw_line) in contents.lines().enumerate() {
            let line = raw_line.split('!').next().unwrap_or("").trim();
            if line.is_empty() || line.starts_with('*') {
                continue;
            }

            let tokens: Vec<&str> = line.split_whitespace().collect();
            if tokens.is_empty() {
                continue;
            }

            let key = tokens[0].to_ascii_uppercase();
            section = match key.as_str() {
                "BOND" | "BONDS" => ParamSection::Bonds,
                "ANGLE" | "ANGLES" | "THETA" => ParamSection::Angles,
                "DIHEDRAL" | "DIHEDRALS" | "PHI" => ParamSection::Dihedrals,
                "IMPROPER" | "IMPROPERS" | "IMPHI" => ParamSection::Impropers,
                "NONBONDED" | "NBONDED" | "NONB" => ParamSection::Nonbonded,
                _ => section,
            };

            if key == "MASS" {
                let atom_type = parse_mass_line(&tokens)
                    .map_err(|e| format!("line {}: {e}", line_number + 1))?;
                ff.atom_types.insert(atom_type.name.clone(), atom_type);
                continue;
            }

            if key == "RESI" || key == "PRES" {
                ff.residue_name = tokens.get(1).map(|v| (*v).to_string());
                continue;
            }

            if key == "ATOM" {
                let mut atom = parse_residue_atom(&tokens)
                    .map_err(|e| format!("line {}: {e}", line_number + 1))?;
                let idx = ff.atoms.len();
                atom.index = idx + 1;
                atom_index_by_name.insert(atom.name.clone(), idx);
                ff.atoms.push(atom);
                continue;
            }

            if key == "BOND" && tokens.len() > 2 && !tokens[1].parse::<f64>().is_ok() {
                parse_topology_pairs(&tokens[1..], &atom_index_by_name, &mut ff.bonds)
                    .map_err(|e| format!("line {}: {e}", line_number + 1))?;
                continue;
            }

            if (key == "ANGL" || key == "ANGLE") && tokens.len() >= 4 {
                parse_topology_triples(&tokens[1..], &atom_index_by_name, &mut ff.angles)
                    .map_err(|e| format!("line {}: {e}", line_number + 1))?;
                continue;
            }

            if (key == "DIHE" || key == "PHI") && tokens.len() >= 5 {
                parse_topology_quads(&tokens[1..], &atom_index_by_name, &mut ff.dihedrals)
                    .map_err(|e| format!("line {}: {e}", line_number + 1))?;
                continue;
            }

            if key == "IMPR" && tokens.len() >= 5 {
                parse_topology_quads(&tokens[1..], &atom_index_by_name, &mut ff.impropers)
                    .map_err(|e| format!("line {}: {e}", line_number + 1))?;
                continue;
            }

            if matches!(section, ParamSection::Bonds)
                && tokens.len() >= 4
                && tokens[2].parse::<f64>().is_ok()
            {
                ff.bond_params.push(parse_bond_param(&tokens)?);
                continue;
            }

            if matches!(section, ParamSection::Angles)
                && tokens.len() >= 5
                && tokens[3].parse::<f64>().is_ok()
            {
                ff.angle_params.push(parse_angle_param(&tokens)?);
                continue;
            }

            if matches!(section, ParamSection::Dihedrals)
                && tokens.len() >= 7
                && tokens[4].parse::<f64>().is_ok()
            {
                ff.dihedral_params.push(parse_dihedral_param(&tokens)?);
                continue;
            }

            if matches!(section, ParamSection::Impropers)
                && tokens.len() >= 7
                && tokens[4].parse::<f64>().is_ok()
            {
                ff.improper_params.push(parse_improper_param(&tokens)?);
                continue;
            }

            if matches!(section, ParamSection::Nonbonded) && tokens.len() >= 4 {
                update_nonbonded_params(&tokens, &mut ff.atom_types)?;
            }
        }

        if ff.atoms.is_empty() {
            return Err(
                "charmm input did not contain residue atom records (ATOM lines)".to_string(),
            );
        }

        Ok(ff)
    }

    pub fn read_file(path: &str) -> Result<Self, String> {
        let contents = fs::read_to_string(path)
            .map_err(|e| format!("failed to read charmm file at '{path}': {e}"))?;
        Self::parse_str(&contents)
    }

    pub fn read_file_new(path: &str) -> Result<Self, String> {
        let contents = match fs::read_to_string(path) {
            Ok(contents) => contents,
            Err(e) => return Err(e.to_string()),
        };
        Self::parse_str(&contents)
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

            let sigma = atom_type
                .sigma
                .ok_or_else(|| format!("missing sigma for atom type '{}'", atom.type_name))?;
            let epsilon = atom_type
                .epsilon
                .ok_or_else(|| format!("missing epsilon for atom type '{}'", atom.type_name))?;

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
                mass: atom.mass.unwrap_or(atom_type.mass),
                energy: 0.0,
                atom_type: idx as f64,
                charge: atom.charge,
            });
        }

        let bonds = self
            .bonds
            .iter()
            .map(|b| {
                let t1 = &self.atoms[b.atom1].type_name;
                let t2 = &self.atoms[b.atom2].type_name;
                let p = self
                    .find_bond_param(t1, t2)
                    .ok_or_else(|| format!("missing bond parameter for {t1}-{t2}"))?;
                Ok(Bond {
                    atom1: b.atom1,
                    atom2: b.atom2,
                    k: p.k,
                    r0: p.r0,
                })
            })
            .collect::<Result<Vec<_>, String>>()?;

        let angles = self
            .angles
            .iter()
            .map(|a| {
                let t1 = &self.atoms[a.atom1].type_name;
                let t2 = &self.atoms[a.atom2].type_name;
                let t3 = &self.atoms[a.atom3].type_name;
                let p = self
                    .find_angle_param(t1, t2, t3)
                    .ok_or_else(|| format!("missing angle parameter for {t1}-{t2}-{t3}"))?;
                Ok(Angle {
                    atom1: a.atom1,
                    atom2: a.atom2,
                    atom3: a.atom3,
                    k: p.k,
                    theta0: p.theta0_deg.to_radians(),
                })
            })
            .collect::<Result<Vec<_>, String>>()?;

        let dihedrals = self
            .dihedrals
            .iter()
            .map(|d| {
                let t1 = &self.atoms[d.atom1].type_name;
                let t2 = &self.atoms[d.atom2].type_name;
                let t3 = &self.atoms[d.atom3].type_name;
                let t4 = &self.atoms[d.atom4].type_name;
                let p = self
                    .find_dihedral_param(t1, t2, t3, t4)
                    .ok_or_else(|| format!("missing dihedral parameter for {t1}-{t2}-{t3}-{t4}"))?;
                Ok(Dihedral {
                    atom1: d.atom1,
                    atom2: d.atom2,
                    atom3: d.atom3,
                    atom4: d.atom4,
                    k: p.k,
                    multiplicity: p.multiplicity,
                    phase: p.phase_deg.to_radians(),
                })
            })
            .collect::<Result<Vec<_>, String>>()?;

        let impropers = self
            .impropers
            .iter()
            .map(|d| {
                let t1 = &self.atoms[d.atom1].type_name;
                let t2 = &self.atoms[d.atom2].type_name;
                let t3 = &self.atoms[d.atom3].type_name;
                let t4 = &self.atoms[d.atom4].type_name;
                let p = self
                    .find_improper_param(t1, t2, t3, t4)
                    .ok_or_else(|| format!("missing improper parameter for {t1}-{t2}-{t3}-{t4}"))?;
                Ok(Improper {
                    atom1: d.atom1,
                    atom2: d.atom2,
                    atom3: d.atom3,
                    atom4: d.atom4,
                    k: p.k,
                    psi0: p.psi0_deg.to_radians(),
                })
            })
            .collect::<Result<Vec<_>, String>>()?;

        Ok(System {
            atoms: particles,
            bonds,
            angles,
            dihedrals,
            impropers,
            ..System::default()
        })
    }

    fn find_bond_param(&self, t1: &str, t2: &str) -> Option<&BondParam> {
        self.bond_params
            .iter()
            .find(|p| (p.t1 == t1 && p.t2 == t2) || (p.t1 == t2 && p.t2 == t1))
    }

    fn find_angle_param(&self, t1: &str, t2: &str, t3: &str) -> Option<&AngleParam> {
        self.angle_params.iter().find(|p| {
            (p.t1 == t1 && p.t2 == t2 && p.t3 == t3) || (p.t1 == t3 && p.t2 == t2 && p.t3 == t1)
        })
    }

    fn find_dihedral_param(
        &self,
        t1: &str,
        t2: &str,
        t3: &str,
        t4: &str,
    ) -> Option<&DihedralParam> {
        self.dihedral_params.iter().find(|p| {
            matches_quad(t1, t2, t3, t4, &p.t1, &p.t2, &p.t3, &p.t4)
                || matches_quad(t4, t3, t2, t1, &p.t1, &p.t2, &p.t3, &p.t4)
        })
    }

    fn find_improper_param(
        &self,
        t1: &str,
        t2: &str,
        t3: &str,
        t4: &str,
    ) -> Option<&ImproperParam> {
        self.improper_params.iter().find(|p| {
            matches_quad(t1, t2, t3, t4, &p.t1, &p.t2, &p.t3, &p.t4)
                || matches_quad(t4, t3, t2, t1, &p.t1, &p.t2, &p.t3, &p.t4)
        })
    }
}

fn parse_mass_line(tokens: &[&str]) -> Result<CharmmAtomType, String> {
    if tokens.len() < 4 {
        return Err("MASS row requires: MASS <id> <type> <mass> [element]".to_string());
    }

    Ok(CharmmAtomType {
        name: tokens[2].to_string(),
        mass: tokens[3]
            .parse::<f64>()
            .map_err(|e| format!("invalid MASS value: {e}"))?,
        element: tokens.get(4).map(|v| (*v).to_string()),
        ..Default::default()
    })
}

fn parse_residue_atom(tokens: &[&str]) -> Result<CharmmAtom, String> {
    if tokens.len() < 4 {
        return Err("ATOM row requires: ATOM <name> <type> <charge> [mass]".to_string());
    }

    Ok(CharmmAtom {
        index: 0,
        name: tokens[1].to_string(),
        type_name: tokens[2].to_string(),
        charge: tokens[3]
            .parse::<f64>()
            .map_err(|e| format!("invalid ATOM charge: {e}"))?,
        mass: tokens.get(4).and_then(|v| v.parse::<f64>().ok()),
    })
}

fn parse_topology_pairs(
    tokens: &[&str],
    atom_index_by_name: &HashMap<String, usize>,
    out: &mut Vec<CharmmBondDef>,
) -> Result<(), String> {
    for pair in tokens.chunks(2) {
        if pair.len() < 2 {
            continue;
        }
        let a1 = atom_index_by_name
            .get(pair[0])
            .ok_or_else(|| format!("unknown atom name '{}'", pair[0]))?;
        let a2 = atom_index_by_name
            .get(pair[1])
            .ok_or_else(|| format!("unknown atom name '{}'", pair[1]))?;
        out.push(CharmmBondDef {
            atom1: *a1,
            atom2: *a2,
        });
    }
    Ok(())
}

fn parse_topology_triples(
    tokens: &[&str],
    atom_index_by_name: &HashMap<String, usize>,
    out: &mut Vec<CharmmAngleDef>,
) -> Result<(), String> {
    for triple in tokens.chunks(3) {
        if triple.len() < 3 {
            continue;
        }
        out.push(CharmmAngleDef {
            atom1: *atom_index_by_name
                .get(triple[0])
                .ok_or_else(|| format!("unknown atom name '{}'", triple[0]))?,
            atom2: *atom_index_by_name
                .get(triple[1])
                .ok_or_else(|| format!("unknown atom name '{}'", triple[1]))?,
            atom3: *atom_index_by_name
                .get(triple[2])
                .ok_or_else(|| format!("unknown atom name '{}'", triple[2]))?,
        });
    }
    Ok(())
}

fn parse_topology_quads(
    tokens: &[&str],
    atom_index_by_name: &HashMap<String, usize>,
    out: &mut Vec<CharmmDihedralDef>,
) -> Result<(), String> {
    for quad in tokens.chunks(4) {
        if quad.len() < 4 {
            continue;
        }
        out.push(CharmmDihedralDef {
            atom1: *atom_index_by_name
                .get(quad[0])
                .ok_or_else(|| format!("unknown atom name '{}'", quad[0]))?,
            atom2: *atom_index_by_name
                .get(quad[1])
                .ok_or_else(|| format!("unknown atom name '{}'", quad[1]))?,
            atom3: *atom_index_by_name
                .get(quad[2])
                .ok_or_else(|| format!("unknown atom name '{}'", quad[2]))?,
            atom4: *atom_index_by_name
                .get(quad[3])
                .ok_or_else(|| format!("unknown atom name '{}'", quad[3]))?,
        });
    }
    Ok(())
}

fn parse_bond_param(tokens: &[&str]) -> Result<BondParam, String> {
    Ok(BondParam {
        t1: tokens[0].to_string(),
        t2: tokens[1].to_string(),
        k: parse_f64(tokens[2], "bond k")?,
        r0: parse_f64(tokens[3], "bond r0")?,
    })
}

fn parse_angle_param(tokens: &[&str]) -> Result<AngleParam, String> {
    Ok(AngleParam {
        t1: tokens[0].to_string(),
        t2: tokens[1].to_string(),
        t3: tokens[2].to_string(),
        k: parse_f64(tokens[3], "angle k")?,
        theta0_deg: parse_f64(tokens[4], "angle theta0")?,
    })
}

fn parse_dihedral_param(tokens: &[&str]) -> Result<DihedralParam, String> {
    Ok(DihedralParam {
        t1: tokens[0].to_string(),
        t2: tokens[1].to_string(),
        t3: tokens[2].to_string(),
        t4: tokens[3].to_string(),
        k: parse_f64(tokens[4], "dihedral k")?,
        multiplicity: parse_usize(tokens[5], "dihedral multiplicity")?,
        phase_deg: parse_f64(tokens[6], "dihedral phase")?,
    })
}

fn parse_improper_param(tokens: &[&str]) -> Result<ImproperParam, String> {
    Ok(ImproperParam {
        t1: tokens[0].to_string(),
        t2: tokens[1].to_string(),
        t3: tokens[2].to_string(),
        t4: tokens[3].to_string(),
        k: parse_f64(tokens[4], "improper k")?,
        psi0_deg: parse_f64(tokens[6], "improper psi0")?,
    })
}

fn update_nonbonded_params(
    tokens: &[&str],
    atom_types: &mut HashMap<String, CharmmAtomType>,
) -> Result<(), String> {
    if tokens[0].eq_ignore_ascii_case("CUTNB")
        || tokens[0].eq_ignore_ascii_case("NBXMOD")
        || tokens[0].eq_ignore_ascii_case("E14FAC")
    {
        return Ok(());
    }

    let type_name = tokens[0];
    let epsilon = tokens[tokens.len() - 2]
        .parse::<f64>()
        .map(|v| v.abs())
        .map_err(|e| format!("invalid nonbonded epsilon: {e}"))?;
    let rmin_half = tokens[tokens.len() - 1]
        .parse::<f64>()
        .map_err(|e| format!("invalid nonbonded Rmin/2: {e}"))?;
    let sigma = rmin_half_to_sigma(rmin_half);

    let atom_type = atom_types
        .entry(type_name.to_string())
        .or_insert_with(|| CharmmAtomType {
            name: type_name.to_string(),
            ..Default::default()
        });
    atom_type.epsilon = Some(epsilon);
    atom_type.sigma = Some(sigma);

    if tokens.len() >= 6 {
        let eps14 = tokens[tokens.len() - 4]
            .parse::<f64>()
            .map(|v| v.abs())
            .ok();
        let rmin14_half = tokens[tokens.len() - 3].parse::<f64>().ok();
        if let (Some(e14), Some(r14)) = (eps14, rmin14_half) {
            atom_type.epsilon_14 = Some(e14);
            atom_type.sigma_14 = Some(rmin_half_to_sigma(r14));
        }
    }

    Ok(())
}

fn rmin_half_to_sigma(rmin_half: f64) -> f64 {
    (2.0 * rmin_half) / 2f64.powf(1.0 / 6.0)
}

fn matches_quad(
    a1: &str,
    a2: &str,
    a3: &str,
    a4: &str,
    p1: &str,
    p2: &str,
    p3: &str,
    p4: &str,
) -> bool {
    matches_type(a1, p1) && matches_type(a2, p2) && matches_type(a3, p3) && matches_type(a4, p4)
}

fn matches_type(value: &str, pattern: &str) -> bool {
    pattern == "X" || pattern == value
}

fn parse_f64(token: &str, what: &str) -> Result<f64, String> {
    token
        .parse::<f64>()
        .map_err(|e| format!("invalid {what}: {e}"))
}

fn parse_usize(token: &str, what: &str) -> Result<usize, String> {
    token
        .parse::<usize>()
        .map_err(|e| format!("invalid {what}: {e}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_charmm_stream_and_builds_system() {
        let input = r#"
* mini charmm stream
*
MASS 1 CT1 12.011 C
MASS 2 HC 1.008 H

RESI MET 0.0
ATOM C1 CT1 -0.18
ATOM H1 HC 0.06
ATOM H2 HC 0.06
ATOM H3 HC 0.06
BOND C1 H1 C1 H2 C1 H3
ANGL H1 C1 H2 H1 C1 H3 H2 C1 H3
DIHE H1 C1 H2 H3
IMPR C1 H1 H2 H3

BONDS
CT1 HC 340.0 1.09

ANGLES
HC CT1 HC 33.0 109.5

DIHEDRALS
X CT1 HC X 0.15 3 0.0

IMPROPER
CT1 HC HC HC 1.1 0 35.0

NONBONDED
CT1 0.0 -0.110 2.00
HC  0.0 -0.020 1.34
"#;

        let ff = CharmmForceField::parse_str(input).expect("charmm parsing should succeed");
        assert_eq!(ff.atoms.len(), 4);
        assert_eq!(ff.bonds.len(), 3);
        assert_eq!(ff.angles.len(), 3);
        assert_eq!(ff.dihedrals.len(), 1);
        assert_eq!(ff.impropers.len(), 1);

        let coords = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];

        let system = ff
            .to_system(&coords)
            .expect("charmm forcefield conversion should succeed");

        assert_eq!(system.atoms.len(), 4);
        assert_eq!(system.bonds.len(), 3);
        assert_eq!(system.angles.len(), 3);
        assert_eq!(system.dihedrals.len(), 1);
        assert_eq!(system.impropers.len(), 1);

        assert!((system.atoms[0].lj_parameters.epsilon - 0.110).abs() < 1e-12);
        let expected_sigma = (2.0 * 2.00) / 2f64.powf(1.0 / 6.0);
        assert!((system.atoms[0].lj_parameters.sigma - expected_sigma).abs() < 1e-12);
    }

    #[test]
    fn errors_when_no_residue_atoms_present() {
        let input = "MASS 1 CT1 12.011 C\n";
        let err = CharmmForceField::parse_str(input).expect_err("should reject missing atoms");
        assert!(err.contains("ATOM"));
    }
}
