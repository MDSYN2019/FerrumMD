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
