use nalgebra::{linalg::SymmetricEigen, DMatrix};

/// Minimal closed-shell Hartree–Fock SCF container.
///
/// All matrices are in an AO basis and use chemist's notation for two-electron
/// integrals: (μν|λσ).

#[derive(Clone, Debug)]
pub struct ScfSystem {
    pub overlap: DMatrix<f64>,
    pub core_hamiltonian: DMatrix<f64>,
    pub eri: Vec<f64>,
    pub n_basis: usize,
    pub n_electrons: usize,
}

#[derive(Clone, Debug)]
pub struct ScfResult {
    pub electronic_energy: f64,
    pub total_energy: f64,
    pub orbital_energies: Vec<f64>,
    pub iterations: usize,
    pub converged: bool,
}

impl ScfSystem {
    fn sort_eigensystem(
        mut eig: SymmetricEigen<f64, nalgebra::Dyn>,
    ) -> SymmetricEigen<f64, nalgebra::Dyn> {
        let n = eig.eigenvalues.len(); // number of eigenvalues/eigenvectors
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| eig.eigenvalues[a].partial_cmp(&eig.eigenvalues[b]).unwrap());

        let mut sorted_vals = eig.eigenvalues.clone();
        let mut sorted_vecs = eig.eigenvectors.clone();

        for (new_col, &old_col) in order.iter().enumerate() {
            sorted_vals[new_col] = eig.eigenvalues[old_col];
            for row in 0..n {
                sorted_vecs[(row, new_col)] = eig.eigenvectors[(row, old_col)];
            }
        }

        eig.eigenvalues = sorted_vals;
        eig.eigenvectors = sorted_vecs;
        eig
    }

    pub fn new(
        overlap: DMatrix<f64>,
        core_hamiltonian: DMatrix<f64>,
        eri: Vec<f64>,
        n_electrons: usize,
    ) -> Self {
        let n_basis = overlap.nrows();
        assert_eq!(overlap.ncols(), n_basis, "overlap must be square");
        assert_eq!(core_hamiltonian.shape(), (n_basis, n_basis));
        assert_eq!(eri.len(), n_basis * n_basis * n_basis * n_basis);

        Self {
            overlap,
            core_hamiltonian,
            eri,
            n_basis,
            n_electrons,
        }
    }

    fn eri_idx(&self, mu: usize, nu: usize, lambda: usize, sigma: usize) -> usize {
        /*
        The ERO indexing is based on a straightforward flattening of the 4D array of integrals
         */

        (((mu * self.n_basis + nu) * self.n_basis + lambda) * self.n_basis) + sigma
    }

    fn eri(&self, mu: usize, nu: usize, lambda: usize, sigma: usize) -> f64 {
        self.eri[self.eri_idx(mu, nu, lambda, sigma)]
    }

    fn orthogonalizer(&self) -> DMatrix<f64> {
        let eig = SymmetricEigen::new(self.overlap.clone());
        let inv_sqrt_vals = DMatrix::from_diagonal(&eig.eigenvalues.map(|x| 1.0 / x.sqrt()));
        &eig.eigenvectors * inv_sqrt_vals * eig.eigenvectors.transpose()
    }

    fn build_fock(&self, density: &DMatrix<f64>) -> DMatrix<f64> {
        let mut fock = self.core_hamiltonian.clone();

        for mu in 0..self.n_basis {
            for nu in 0..self.n_basis {
                let mut g_mu_nu = 0.0;
                for lambda in 0..self.n_basis {
                    for sigma in 0..self.n_basis {
                        let coulomb = self.eri(mu, nu, lambda, sigma);
                        let exchange = self.eri(mu, lambda, nu, sigma);
                        g_mu_nu += density[(lambda, sigma)] * (coulomb - 0.5 * exchange);
                    }
                }
                fock[(mu, nu)] += g_mu_nu;
            }
        }

        fock
    }

    fn build_density(&self, coeff: &DMatrix<f64>) -> DMatrix<f64> {
        let n_occ = self.n_electrons / 2;
        let mut density = DMatrix::zeros(self.n_basis, self.n_basis);

        for mu in 0..self.n_basis {
            for nu in 0..self.n_basis {
                let mut value = 0.0;
                for m in 0..n_occ {
                    value += coeff[(mu, m)] * coeff[(nu, m)];
                }
                density[(mu, nu)] = 2.0 * value;
            }
        }

        density
    }

    fn electronic_energy(&self, density: &DMatrix<f64>, fock: &DMatrix<f64>) -> f64 {
        let mut energy = 0.0;
        for mu in 0..self.n_basis {
            for nu in 0..self.n_basis {
                energy += density[(mu, nu)] * (self.core_hamiltonian[(mu, nu)] + fock[(mu, nu)]);
            }
        }
        0.5 * energy
    }

    pub fn run_scf(
        &self,
        nuclear_repulsion: f64,
        max_iter: usize,
        energy_tol: f64,
        density_tol: f64,
    ) -> ScfResult {
        let x = self.orthogonalizer();

        let mut density = DMatrix::zeros(self.n_basis, self.n_basis);
        let mut old_energy = f64::INFINITY;
        let mut orbital_energies = vec![0.0; self.n_basis];

        for iter in 1..=max_iter {
            let fock = self.build_fock(&density);
            let fock_ortho = x.transpose() * &fock * &x;
            let eig = Self::sort_eigensystem(SymmetricEigen::new(fock_ortho));
            let coeff = &x * eig.eigenvectors;
            let new_density = self.build_density(&coeff);

            let energy = self.electronic_energy(&new_density, &fock);
            let d_energy = (energy - old_energy).abs();
            let d_density = (&new_density - &density).norm();

            orbital_energies.clone_from_slice(eig.eigenvalues.as_slice());

            if d_energy < energy_tol && d_density < density_tol {
                return ScfResult {
                    electronic_energy: energy,
                    total_energy: energy + nuclear_repulsion,
                    orbital_energies,
                    iterations: iter,
                    converged: true,
                };
            }

            density = new_density;
            old_energy = energy;
        }

        let fock = self.build_fock(&density);
        let energy = self.electronic_energy(&density, &fock);

        ScfResult {
            electronic_energy: energy,
            total_energy: energy + nuclear_repulsion,
            orbital_energies,
            iterations: max_iter,
            converged: false,
        }
    }
}

/// A fragment definition for fragment molecular orbital (FMO) calculations.
///
/// `basis_indices` selects the AO rows/columns owned by the fragment in the
/// parent `ScfSystem`. `n_electrons` must currently be even because the
/// underlying SCF implementation is restricted to closed-shell RHF.
#[derive(Clone, Debug)]
pub struct FmoFragment {
    pub name: String,
    pub basis_indices: Vec<usize>,
    pub n_electrons: usize,
    pub nuclear_repulsion: f64,
    pub net_charge: f64,
    pub center: [f64; 3],
}

impl FmoFragment {
    pub fn new(
        name: impl Into<String>,
        basis_indices: Vec<usize>,
        n_electrons: usize,
        nuclear_repulsion: f64,
    ) -> Self {
        Self {
            name: name.into(),
            basis_indices,
            n_electrons,
            nuclear_repulsion,
            net_charge: 0.0,
            center: [0.0; 3],
        }
    }

    pub fn with_pieda_site(mut self, net_charge: f64, center: [f64; 3]) -> Self {
        self.net_charge = net_charge;
        self.center = center;
        self
    }
}

#[derive(Clone, Debug)]
pub struct FmoPairConfig {
    pub i: usize,
    pub j: usize,
    pub nuclear_repulsion: f64,
}

#[derive(Clone, Debug)]
pub struct FmoMonomerResult {
    pub fragment: String,
    pub energy: f64,
    pub scf: ScfResult,
}

#[derive(Clone, Debug)]
pub struct PiedaPairResult {
    pub fragment_i: String,
    pub fragment_j: String,
    pub dimer_energy: f64,
    pub interaction_energy: f64,
    pub electrostatic: f64,
    pub exchange_repulsion: f64,
    pub charge_transfer: f64,
    pub dispersion: f64,
}

#[derive(Clone, Debug)]
pub struct FmoResult {
    pub total_energy: f64,
    pub monomers: Vec<FmoMonomerResult>,
    pub pairs: Vec<PiedaPairResult>,
    pub converged: bool,
}

/// FMO2 driver backed by the closed-shell Hartree-Fock SCF implementation.
///
/// This implementation computes embedded-data-free gas-phase FMO2 energies:
/// E(FMO2) = Σ_I E_I + Σ_{I<J}(E_IJ - E_I - E_J).  The PIEDA fields are a
/// pragmatic decomposition suitable for diagnostics: electrostatics is computed
/// from fragment point charges and centers, exchange-repulsion is the remaining
/// HF interaction after optional charge-transfer/dispersion placeholders.
pub struct FmoSystem {
    pub parent: ScfSystem,
    pub fragments: Vec<FmoFragment>,
    pub pair_configs: Vec<FmoPairConfig>,
}

impl FmoSystem {
    pub fn new(parent: ScfSystem, fragments: Vec<FmoFragment>) -> Self {
        let mut pair_configs = Vec::new();
        for i in 0..fragments.len() {
            for j in (i + 1)..fragments.len() {
                pair_configs.push(FmoPairConfig {
                    i,
                    j,
                    nuclear_repulsion: 0.0,
                });
            }
        }
        Self {
            parent,
            fragments,
            pair_configs,
        }
    }

    pub fn with_pair_nuclear_repulsion(
        mut self,
        i: usize,
        j: usize,
        nuclear_repulsion: f64,
    ) -> Self {
        if let Some(pair) = self
            .pair_configs
            .iter_mut()
            .find(|p| (p.i == i && p.j == j) || (p.i == j && p.j == i))
        {
            pair.nuclear_repulsion = nuclear_repulsion;
        }
        self
    }

    fn sub_system(&self, basis_indices: &[usize], n_electrons: usize) -> ScfSystem {
        let n = basis_indices.len();
        let overlap = DMatrix::from_fn(n, n, |a, b| {
            self.parent.overlap[(basis_indices[a], basis_indices[b])]
        });
        let core_hamiltonian = DMatrix::from_fn(n, n, |a, b| {
            self.parent.core_hamiltonian[(basis_indices[a], basis_indices[b])]
        });
        let mut eri = vec![0.0; n * n * n * n];
        for mu in 0..n {
            for nu in 0..n {
                for lambda in 0..n {
                    for sigma in 0..n {
                        let sub_idx = (((mu * n + nu) * n + lambda) * n) + sigma;
                        eri[sub_idx] = self.parent.eri(
                            basis_indices[mu],
                            basis_indices[nu],
                            basis_indices[lambda],
                            basis_indices[sigma],
                        );
                    }
                }
            }
        }
        ScfSystem::new(overlap, core_hamiltonian, eri, n_electrons)
    }

    fn pair_basis(&self, i: usize, j: usize) -> Vec<usize> {
        let mut basis = self.fragments[i].basis_indices.clone();
        basis.extend_from_slice(&self.fragments[j].basis_indices);
        basis
    }

    fn point_charge_electrostatic(&self, i: usize, j: usize) -> f64 {
        let fi = &self.fragments[i];
        let fj = &self.fragments[j];
        let dx = fi.center[0] - fj.center[0];
        let dy = fi.center[1] - fj.center[1];
        let dz = fi.center[2] - fj.center[2];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        if r <= f64::EPSILON {
            0.0
        } else {
            fi.net_charge * fj.net_charge / r
        }
    }

    pub fn run_fmo2(&self, max_iter: usize, energy_tol: f64, density_tol: f64) -> FmoResult {
        let mut monomers = Vec::with_capacity(self.fragments.len());
        let mut converged = true;
        for fragment in &self.fragments {
            let system = self.sub_system(&fragment.basis_indices, fragment.n_electrons);
            let scf = system.run_scf(
                fragment.nuclear_repulsion,
                max_iter,
                energy_tol,
                density_tol,
            );
            converged &= scf.converged;
            monomers.push(FmoMonomerResult {
                fragment: fragment.name.clone(),
                energy: scf.total_energy,
                scf,
            });
        }

        let mut total_energy: f64 = monomers.iter().map(|m| m.energy).sum();
        let mut pairs = Vec::with_capacity(self.pair_configs.len());
        for pair in &self.pair_configs {
            let basis = self.pair_basis(pair.i, pair.j);
            let electrons = self.fragments[pair.i].n_electrons + self.fragments[pair.j].n_electrons;
            let scf = self.sub_system(&basis, electrons).run_scf(
                pair.nuclear_repulsion,
                max_iter,
                energy_tol,
                density_tol,
            );
            converged &= scf.converged;
            let interaction_energy =
                scf.total_energy - monomers[pair.i].energy - monomers[pair.j].energy;
            total_energy += interaction_energy;
            let electrostatic = self.point_charge_electrostatic(pair.i, pair.j);
            let charge_transfer = 0.0;
            let dispersion = 0.0;
            let exchange_repulsion =
                interaction_energy - electrostatic - charge_transfer - dispersion;
            pairs.push(PiedaPairResult {
                fragment_i: self.fragments[pair.i].name.clone(),
                fragment_j: self.fragments[pair.j].name.clone(),
                dimer_energy: scf.total_energy,
                interaction_energy,
                electrostatic,
                exchange_repulsion,
                charge_transfer,
                dispersion,
            });
        }

        FmoResult {
            total_energy,
            monomers,
            pairs,
            converged,
        }
    }
}

use std::cmp::Ordering;
use std::collections::HashMap;
use std::error::Error;
use std::fmt;
use std::io::BufRead;

/// One row from the FMO-SCOP-29Jun2022 IFIE/PIEDA TSV files described by
/// Takaya et al. Scientific Data 2024.
#[derive(Clone, Debug, PartialEq)]
pub struct FmoScopRecord {
    pub pdb_id: String,
    pub i_chain: String,
    pub j_chain: String,
    pub i_pair: i32,
    pub j_pair: i32,
    pub i_name: String,
    pub j_name: String,
    pub ifie: f64,
    pub electrostatic: f64,
    pub exchange_repulsion: f64,
    pub charge_transfer_mix: f64,
    pub dispersion: f64,
    pub distance: f64,
    pub approx: bool,
    pub i_charge: i32,
    pub j_charge: i32,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FmoScopParseError {
    pub line: usize,
    pub message: String,
}

impl fmt::Display for FmoScopParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "FMO-SCOP parse error on line {}: {}",
            self.line, self.message
        )
    }
}

impl Error for FmoScopParseError {}

#[derive(Copy, Clone, Debug, Eq, Hash, PartialEq)]
pub enum FragmentChargeClass {
    NeutralNeutral,
    AttractiveCharged,
    PositiveNeutral,
    NegativeNeutral,
    PositivePositive,
    NegativeNegative,
    Other,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PiedaSummary {
    pub count: usize,
    pub median_ifie: f64,
    pub median_electrostatic: f64,
    pub median_exchange_repulsion: f64,
    pub median_charge_transfer_mix: f64,
    pub median_dispersion: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ResiduePairPiedaSummary {
    pub residue_a: String,
    pub residue_b: String,
    pub summary: PiedaSummary,
}

impl FmoScopRecord {
    pub fn parse_tsv_line(line_number: usize, line: &str) -> Result<Self, FmoScopParseError> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 16 {
            return Err(FmoScopParseError {
                line: line_number,
                message: format!("expected 16 tab-separated fields, found {}", fields.len()),
            });
        }

        Ok(Self {
            pdb_id: fields[0].to_string(),
            i_chain: fields[1].to_string(),
            j_chain: fields[2].to_string(),
            i_pair: parse_field(line_number, "Ipair", fields[3])?,
            j_pair: parse_field(line_number, "Jpair", fields[4])?,
            i_name: fields[5].to_string(),
            j_name: fields[6].to_string(),
            ifie: parse_field(line_number, "IFIE", fields[7])?,
            electrostatic: parse_field(line_number, "ES", fields[8])?,
            exchange_repulsion: parse_field(line_number, "EX", fields[9])?,
            charge_transfer_mix: parse_field(line_number, "CT", fields[10])?,
            dispersion: parse_field(line_number, "DI_MP2", fields[11])?,
            distance: parse_field(line_number, "Dist", fields[12])?,
            approx: match fields[13] {
                "T" | "t" | "true" | "TRUE" | "1" => true,
                "F" | "f" | "false" | "FALSE" | "0" => false,
                value => {
                    return Err(FmoScopParseError {
                        line: line_number,
                        message: format!("invalid approx value '{value}'"),
                    })
                }
            },
            i_charge: parse_field(line_number, "Ielect", fields[14])?,
            j_charge: parse_field(line_number, "Jelect", fields[15])?,
        })
    }

    pub fn charge_class(&self) -> FragmentChargeClass {
        match (self.i_charge.signum(), self.j_charge.signum()) {
            (0, 0) => FragmentChargeClass::NeutralNeutral,
            (1, -1) | (-1, 1) => FragmentChargeClass::AttractiveCharged,
            (1, 0) | (0, 1) => FragmentChargeClass::PositiveNeutral,
            (-1, 0) | (0, -1) => FragmentChargeClass::NegativeNeutral,
            (1, 1) => FragmentChargeClass::PositivePositive,
            (-1, -1) => FragmentChargeClass::NegativeNegative,
            _ => FragmentChargeClass::Other,
        }
    }
}

fn parse_field<T>(line: usize, field: &str, value: &str) -> Result<T, FmoScopParseError>
where
    T: std::str::FromStr,
    T::Err: fmt::Display,
{
    value.parse::<T>().map_err(|err| FmoScopParseError {
        line,
        message: format!("invalid {field} value '{value}': {err}"),
    })
}

pub fn read_fmo_scop_tsv<R: BufRead>(reader: R) -> Result<Vec<FmoScopRecord>, FmoScopParseError> {
    let mut records = Vec::new();
    for (idx, line) in reader.lines().enumerate() {
        let line_number = idx + 1;
        let line = line.map_err(|err| FmoScopParseError {
            line: line_number,
            message: err.to_string(),
        })?;
        let trimmed = line.trim_end();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        if records.is_empty() && looks_like_header(trimmed) {
            continue;
        }
        records.push(FmoScopRecord::parse_tsv_line(line_number, trimmed)?);
    }
    Ok(records)
}

fn looks_like_header(line: &str) -> bool {
    let first = line.split('\t').next().unwrap_or_default();
    first.eq_ignore_ascii_case("pdbid") || first.eq_ignore_ascii_case("pdb_id")
}

pub fn summarize_by_charge_class(
    records: &[FmoScopRecord],
) -> HashMap<FragmentChargeClass, PiedaSummary> {
    let mut grouped: HashMap<FragmentChargeClass, Vec<&FmoScopRecord>> = HashMap::new();
    for record in records.iter().filter(|record| !record.approx) {
        grouped
            .entry(record.charge_class())
            .or_default()
            .push(record);
    }
    grouped
        .into_iter()
        .map(|(class, values)| (class, summarize_records(&values)))
        .collect()
}

pub fn summarize_by_residue_pair(records: &[FmoScopRecord]) -> Vec<ResiduePairPiedaSummary> {
    let mut grouped: HashMap<(String, String), Vec<&FmoScopRecord>> = HashMap::new();
    for record in records.iter().filter(|record| !record.approx) {
        let mut pair = [record.i_name.clone(), record.j_name.clone()];
        pair.sort();
        grouped
            .entry((pair[0].clone(), pair[1].clone()))
            .or_default()
            .push(record);
    }
    let mut summaries: Vec<_> = grouped
        .into_iter()
        .map(|((residue_a, residue_b), values)| ResiduePairPiedaSummary {
            residue_a,
            residue_b,
            summary: summarize_records(&values),
        })
        .collect();
    summaries.sort_by(|a, b| match a.residue_a.cmp(&b.residue_a) {
        Ordering::Equal => a.residue_b.cmp(&b.residue_b),
        other => other,
    });
    summaries
}

fn summarize_records(records: &[&FmoScopRecord]) -> PiedaSummary {
    PiedaSummary {
        count: records.len(),
        median_ifie: median(records.iter().map(|record| record.ifie).collect()),
        median_electrostatic: median(records.iter().map(|record| record.electrostatic).collect()),
        median_exchange_repulsion: median(
            records
                .iter()
                .map(|record| record.exchange_repulsion)
                .collect(),
        ),
        median_charge_transfer_mix: median(
            records
                .iter()
                .map(|record| record.charge_transfer_mix)
                .collect(),
        ),
        median_dispersion: median(records.iter().map(|record| record.dispersion).collect()),
    }
}

fn median(mut values: Vec<f64>) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    values.sort_by(|a, b| a.total_cmp(b));
    let mid = values.len() / 2;
    if values.len() % 2 == 0 {
        (values[mid - 1] + values[mid]) / 2.0
    } else {
        values[mid]
    }
}

/// Shared output object for SCF-like quantum templates.
pub struct TemplateScfResult {
    pub electronic_energy: f64,
    pub total_energy: f64,
    pub iterations: usize,
    pub converged: bool,
    pub details: String,
}

/// A lightweight Hartree–Fock template that reuses the existing `ScfSystem`
/// implementation for integrals and SCF controls.
pub struct HartreeFockTemplate {
    pub system: ScfSystem,
}

impl HartreeFockTemplate {
    pub fn new(system: ScfSystem) -> Self {
        Self { system }
    }

    pub fn run_template(
        &self,
        nuclear_repulsion: f64,
        max_iter: usize,
        energy_tol: f64,
        density_tol: f64,
    ) -> TemplateScfResult {
        let result = self
            .system
            .run_scf(nuclear_repulsion, max_iter, energy_tol, density_tol);

        TemplateScfResult {
            electronic_energy: result.electronic_energy,
            total_energy: result.total_energy,
            iterations: result.iterations,
            converged: result.converged,
            details: "HF template: closed-shell SCF using Coulomb + exchange terms".to_string(),
        }
    }
}

/// A starter Kohn–Sham DFT template that mirrors the HF data flow and keeps
/// explicit placeholders for exchange-correlation pieces.
pub struct DftTemplate {
    pub system: ScfSystem,
    pub xc_functional: String,
    pub grid_points: usize,
}

impl DftTemplate {
    pub fn new(system: ScfSystem, xc_functional: impl Into<String>, grid_points: usize) -> Self {
        Self {
            system,
            xc_functional: xc_functional.into(),
            grid_points,
        }
    }

    /// Runs the DFT template.
    ///
    /// TODO:
    /// - Add numerical integration over atom-centered grids.
    /// - Build V_xc and E_xc from `xc_functional`.
    /// - Replace the HF exchange term with KS exchange-correlation treatment.
    pub fn run_template(
        &self,
        nuclear_repulsion: f64,
        max_iter: usize,
        energy_tol: f64,
        density_tol: f64,
    ) -> TemplateScfResult {
        let hf_like_result =
            self.system
                .run_scf(nuclear_repulsion, max_iter, energy_tol, density_tol);

        TemplateScfResult {
            electronic_energy: hf_like_result.electronic_energy,
            total_energy: hf_like_result.total_energy,
            iterations: hf_like_result.iterations,
            converged: false,
            details: format!(
                "DFT template placeholder with XC='{}' and {} grid points; replace HF exchange with V_xc",
                self.xc_functional, self.grid_points
            ),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hartree_fock_single_basis_closed_shell() {
        let overlap = DMatrix::from_row_slice(1, 1, &[1.0]);
        let core_h = DMatrix::from_row_slice(1, 1, &[-1.0]);
        let eri = vec![0.7]; // (00|00)

        let hf = ScfSystem::new(overlap, core_h, eri, 2);
        let result = hf.run_scf(0.0, 50, 1e-12, 1e-12);

        assert!(result.converged);
        assert!(result.iterations <= 3);
        assert!((result.electronic_energy + 1.3).abs() < 1e-10);
        assert!((result.total_energy + 1.3).abs() < 1e-10);
        assert!((result.orbital_energies[0] + 0.3).abs() < 1e-10);
    }

    #[test]
    fn hartree_fock_two_basis_non_interacting_case() {
        // Non-interacting closed-shell test: ERIs are all zero, so SCF should
        // reduce to diagonalizing the core Hamiltonian in an orthonormal basis.
        let overlap = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);
        let core_h = DMatrix::from_row_slice(2, 2, &[-1.0, -0.2, -0.2, -0.5]);
        let eri = vec![0.0; 16];

        let hf = ScfSystem::new(overlap, core_h, eri, 2);
        let result = hf.run_scf(0.0, 50, 1e-12, 1e-12);

        assert!(result.converged);
        assert!(result.iterations <= 3);

        // Lowest eigenvalue of core_h is ~ -1.0701562118716426 for this matrix.
        // With a closed shell (2 electrons), E_elec = 2 * e_occ.
        let expected_electronic = -2.140312423743285;
        assert!((result.electronic_energy - expected_electronic).abs() < 1e-10);
        assert!((result.total_energy - expected_electronic).abs() < 1e-10);
        assert!(result.orbital_energies[0] < result.orbital_energies[1]);
    }

    #[test]
    fn fmo2_and_pieda_two_fragment_non_interacting_case() {
        let overlap = DMatrix::identity(2, 2);
        let core_h = DMatrix::from_diagonal(&nalgebra::DVector::from_vec(vec![-1.0, -0.5]));
        let eri = vec![0.0; 16];
        let parent = ScfSystem::new(overlap, core_h, eri, 4);
        let fragments = vec![
            FmoFragment::new("A", vec![0], 2, 0.0).with_pieda_site(1.0, [0.0, 0.0, 0.0]),
            FmoFragment::new("B", vec![1], 2, 0.0).with_pieda_site(-1.0, [2.0, 0.0, 0.0]),
        ];

        let fmo = FmoSystem::new(parent, fragments);
        let result = fmo.run_fmo2(50, 1e-12, 1e-12);

        assert!(result.converged);
        assert_eq!(result.monomers.len(), 2);
        assert_eq!(result.pairs.len(), 1);
        assert!((result.monomers[0].energy + 2.0).abs() < 1e-10);
        assert!((result.monomers[1].energy + 1.0).abs() < 1e-10);
        assert!((result.total_energy + 3.0).abs() < 1e-10);
        assert!((result.pairs[0].interaction_energy).abs() < 1e-10);
        assert!((result.pairs[0].electrostatic + 0.5).abs() < 1e-10);
        assert!((result.pairs[0].exchange_repulsion - 0.5).abs() < 1e-10);
    }

    #[test]
    fn fmo_scop_tsv_parser_and_charge_summaries() {
        let data = "pdbid\tIchain\tJchain\tIpair\tJpair\tIname\tJname\tIFIE\tES\tEX\tCT\tDI_MP2\tDist\tapprox\tIelect\tJelect\n\
1abc\tA\tA\t1\t2\tASP\tLYS\t-10.0\t-12.0\t3.0\t-0.5\t-0.5\t2.0\tF\t-1\t1\n\
1abc\tA\tA\t3\t4\tALA\tGLY\t-2.0\t-1.0\t0.5\t-0.2\t-0.1\t4.0\tF\t0\t0\n\
1abc\tA\tA\t5\t6\tALA\tGLY\t-100.0\t-100.0\t0.0\t0.0\t0.0\t9.0\tT\t0\t0\n";

        let records = read_fmo_scop_tsv(std::io::Cursor::new(data)).unwrap();
        assert_eq!(records.len(), 3);
        assert_eq!(
            records[0].charge_class(),
            FragmentChargeClass::AttractiveCharged
        );
        assert!(records[2].approx);

        let charge_summaries = summarize_by_charge_class(&records);
        let neutral = charge_summaries
            .get(&FragmentChargeClass::NeutralNeutral)
            .expect("neutral-neutral summary");
        assert_eq!(neutral.count, 1);
        assert_eq!(neutral.median_ifie, -2.0);

        let attractive = charge_summaries
            .get(&FragmentChargeClass::AttractiveCharged)
            .expect("attractive charged summary");
        assert_eq!(attractive.count, 1);
        assert_eq!(attractive.median_electrostatic, -12.0);
    }

    #[test]
    fn fmo_scop_residue_pair_summaries_are_order_independent() {
        let data = "1abc\tA\tA\t1\t2\tGLY\tALA\t-1.0\t-0.5\t0.2\t-0.1\t-0.1\t3.0\tF\t0\t0\n\
1abc\tA\tA\t3\t4\tALA\tGLY\t-3.0\t-1.5\t0.6\t-0.3\t-0.3\t3.5\tF\t0\t0\n";

        let records = read_fmo_scop_tsv(std::io::Cursor::new(data)).unwrap();
        let summaries = summarize_by_residue_pair(&records);
        assert_eq!(summaries.len(), 1);
        assert_eq!(summaries[0].residue_a, "ALA");
        assert_eq!(summaries[0].residue_b, "GLY");
        assert_eq!(summaries[0].summary.count, 2);
        assert_eq!(summaries[0].summary.median_ifie, -2.0);
        assert_eq!(summaries[0].summary.median_exchange_repulsion, 0.4);
    }

    #[test]
    fn templates_for_hf_and_dft_are_constructible() {
        let overlap = DMatrix::from_row_slice(1, 1, &[1.0]);
        let core_h = DMatrix::from_row_slice(1, 1, &[-1.0]);
        let eri = vec![0.7];

        let hf_system = ScfSystem::new(overlap.clone(), core_h.clone(), eri.clone(), 2);
        let hf_template = HartreeFockTemplate::new(hf_system);
        let hf_result = hf_template.run_template(0.0, 25, 1e-12, 1e-12);
        assert!(hf_result.converged);

        let dft_system = ScfSystem::new(overlap, core_h, eri, 2);
        let dft_template = DftTemplate::new(dft_system, "PBE", 302);
        let dft_result = dft_template.run_template(0.0, 25, 1e-12, 1e-12);
        assert!(!dft_result.converged);
        assert!(dft_result.details.contains("PBE"));
    }
}
