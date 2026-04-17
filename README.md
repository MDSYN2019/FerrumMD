
![image info](FerrumMD.png)


A high-performance **Molecular Dynamics (MD) engine written in Rust**.

Implements Lennard-Jones (LJ) interactions, bonded forces, velocity-Verlet time integration, periodic boundary conditions, thermostats, and molecular systems such as H₂.

---

# FerrumMD
![CI](https://github.com/MDSYN2019/rust_md_work/actions/workflows/ci.yml/badge.svg)

## 🔥 Features

### ✅ Core MD Functionality
- Velocity Verlet integrator  
- Periodic Boundary Conditions (PBC)  
- Minimum Image Convention  
- Site–site Lennard-Jones interactions  
- Bonded interactions via harmonic springs  
- Support for both:
  - **Particle collections** (`InitOutput::Particles`)
  - **Molecular systems** (`InitOutput::Systems`)  

### 🌡 Thermostat Algorithms
- Maxwell–Boltzmann initial velocity sampling  
- Berendsen thermostat (velocity rescaling)  
- Nose-Hoover thermostat (extended-systems coupling)
- Nose-Hoover isotropic barostat (volume/position scaling)
- NVE and pseudo-NVT control  
- Temperature calculation from kinetic energy  

### 🧬 Molecular Support
- Construction of small molecules (e.g., H₂)  
- Support for multiple molecules via system cloning  
- Bonded forces with equilibrium distances and spring constants  

### 🧰 Utilities
- Energy reporting (kinetic, potential, total)  
- Force calculations (bonded + nonbonded LJ)  
- pbc wrapping  
- Configurable time-step, LJ parameters, masses, and box sizes  
- PDB and GRO coordinate readers (`molecule::io::{read_pdb, read_gro}`)  
- LAMMPS data/dump readers (`molecule::io::{read_lammps_data, read_lammps_dump}`)  
- Martini `.itp` force-field reader + converter (`molecule::martini::MartiniForceField`)  
- CHARMM-style topology/parameter reader + converter (`molecule::charmm::CharmmForceField`)  

---

## 📏 Unit conventions (explicit)

FerrumMD uses a single explicit MD unit system internally:

- **Coordinates / box lengths:** `nm`
- **Velocities:** `nm/ps`
- **Forces:** `kJ/(mol·nm)`
- **Charges:** elementary charge `e`
- **Masses:** atomic mass unit `amu`
- **Timestep (`dt`):** `ps`
- **Energies (kinetic, LJ, Coulomb, total):** `kJ/mol`
- **Coulomb prefactor:** `138.935457644 kJ mol^-1 nm e^-2`

GRO I/O is now unit-preserving (nm in file and nm in memory), so there is no hidden nm↔Å conversion path.

---

## 🚀 Getting Started

### Install Rust
You’ll need a stable Rust toolchain:

```bash
rustup update


## ⚡ MPI Parallel NVE Example

An MPI-enabled NVE integration path is available behind the `mpi` feature flag.
It parallelizes Lennard-Jones force/energy accumulation across ranks and uses collective reductions to build global forces.

```bash
cargo run --features mpi
mpirun -n 4 cargo run --features mpi
```

The MPI code path is intended as a parallel-programming example (`run_md_nve_mpi` / `run_md_nve_particles_mpi`).


### Make targets (serial + MPI)

```bash
make run
make run-mpi NP=4
```

## 🧪 Point-particle water-box style trajectory output (GRO + XTC)

This repository now includes a runnable example binary that creates a simple Lennard-Jones point-particle fluid and writes outputs that can be opened in VMD:

```bash
cargo run --bin water_box
```

Generated files:
- `water_box.gro` (final frame structure)
- `water_box.xtc` (trajectory)

### Open in VMD
1. `vmd water_box.gro`
2. In the VMD GUI: **File → Load Data Into Molecule...**
3. Select `water_box.xtc` and load.

> Notes:
> - This is a coarse-grained, point-particle fluid setup (water-like in mass/density intent, not explicit 3-site/4-site water geometry).
> - Coordinates are written in GRO/XTC-compatible units (nm in files).

## 🖥️ Rust VMD-style trajectory viewer (workbench starter)

A new interactive viewer binary is included as a Rust-first starting point for trajectory visualization workflows.
It loads a `GRO` structure and optionally an `XTC` trajectory, then renders atoms in a 3D viewport with playback controls.

Run it with:

```bash
cargo run --features viewer --bin vmd_viewer -- water_box.gro water_box.xtc
```

Controls:
- `Space`: play/pause
- `←` / `→`: step one frame
- `↑` / `↓`: playback FPS
- Left mouse drag: orbit camera
- Mouse wheel: zoom

This is intentionally a scaffold so you can extend it in Rust (atom coloring, bonds, selection, measure tools, rendering styles, etc.).

---

## 🐍 Python interface (buildable scaffold)

A minimal Python extension interface is available behind the `python` feature.
It exposes:

- `PyMdEngine(sigma, epsilon)` class
- `PyMdEngine.force_at_distance(r)`
- `lj_force_scalar(r, sigma, epsilon)`
- `python_api_version()`

### Build with maturin

```bash
pip install maturin
maturin develop --features python
```

Then in Python:

```python
import sang_md_py

engine = sang_md_py.PyMdEngine(1.0, 1.0)
print(engine.force_at_distance(1.2))
print(sang_md_py.lj_force_scalar(1.2, 1.0, 1.0))
```

This is intended as a starting point you can expand with trajectory stepping, system builders, and observables.

## 🧪 Martini coarse-grained water-box NVT example

A dedicated Martini-style coarse-grained water box example is also available. It runs a single-bead solvent in an NVT-like setup using velocity-Verlet integration with a Berendsen thermostat and writes GRO/XTC outputs:

```bash
cargo run --bin martini_water_box
```

Generated files:
- `martini_water_box.gro`
- `martini_water_box.xtc`

This is intended as a lightweight CG solvent demo that you can visualize in VMD with the same loading flow used for `water_box.xtc`.

---

## 🗺️ TIP3P roadmap to full simulation fidelity

This project includes a TIP3P-oriented water-box path, and the following milestones track the remaining work needed to reach full, production-grade TIP3P simulation fidelity.

### M1 — Rigid water support (highest priority)
**Goal:** Ensure physically correct rigid water geometry and stable integration behavior.

**Scope**
- Parse and use `[ settles ]` water constraints from ITP.
- Implement rigid-constraint enforcement (prefer SETTLE for 3-site water; optional SHAKE/RATTLE fallback).
- Keep flexible and rigid water modes explicit in configuration.

**Acceptance criteria**
- O–H and H–H distances remain within tight tolerance during long runs.
- Stable with larger practical time steps than unconstrained flexible mode.

### M2 — Topology-correct nonbonded handling
**Goal:** Correctly apply molecular exclusions and pair handling.

**Scope**
- Parse and store `[ exclusions ]` (and optionally `[ pairs ]`) topology sections.
- Extend the system representation with exclusion/pair metadata.
- Apply exclusion/pair logic in LJ and electrostatic kernels.

**Acceptance criteria**
- Unit tests confirm excluded pairs do not contribute to nonbonded terms.
- Energy decomposition reflects intended topology behavior.

### M3 — Units and electrostatics audit
**Goal:** Make electrostatics and parameter units unambiguous and physically consistent.

**Scope**
- Replace reduced-unit shortcuts with explicit, documented unit-aware handling.
- Document units for coordinates, charges, mass, timestep, and energy.
- Add regression checks for analytic charge-pair and small water-cluster cases.

**Acceptance criteria**
- Pairwise reference checks pass.
- No hidden or implicit unit conversion paths remain.

### M4 — Integrator and ensemble hygiene
**Goal:** Improve simulation protocol clarity and reproducibility.

**Scope**
- Clarify NVE/NVT naming and behavior in simulation entry points.
- Add explicit equilibration vs production phases in TIP3P examples.
- Remove hidden thermostat assumptions in core loops.

**Acceptance criteria**
- Reproducible runs with a fixed RNG seed option.
- Clear runtime configuration for thermostat/barostat behavior.

### M5 — Performance optimization for realistic box sizes
**Goal:** Achieve practical throughput for larger water systems.

**Scope**
- Profile and optimize electrostatics (real + reciprocal) and neighbor updates.
- Reduce avoidable temporary allocations and redundant data transforms.
- Add benchmark cases for representative water-box sizes (e.g., 216/512/1024 waters).

**Acceptance criteria**
- Throughput metrics are documented and tracked.
- TIP3P examples complete in practical wall-clock time at target sizes.

### M6 — Scientific validation suite
**Goal:** Validate output quality against expected TIP3P behavior.

**Scope**
- Add observables: density, O–O RDF, diffusion estimate, and NVE energy drift.
- Compare against known TIP3P reference ranges under matching conditions.
- Publish a validation report (markdown + plots).

**Acceptance criteria**
- Validation report is versioned in the repository.
- CI smoke checks cover a reduced-size validation subset.

### Suggested PR sequence
1. ITP parser extensions (`settles`, `exclusions`, optional `pairs`) + data model updates.
2. SETTLE/constraint implementation + focused tests.
3. Exclusion/pair logic integrated into nonbonded calculations.
4. Unit and electrostatics normalization with regression tests.
5. TIP3P protocol cleanup (equilibration/production config).
6. Performance tuning and benchmark documentation.
7. Validation metrics/report integration.

---

## 🧪 Computational Physics for Drug Discovery: high-impact features to add

To support an end-to-end computational chemistry pipeline (hit finding → lead optimization), FerrumMD can be extended with the following capabilities.

### 1) Free-energy platform (ABFE/RBFE)
- Alchemical transformation engine (`λ` windows, soft-core LJ/Coulomb, decoupling/annihilation modes).
- MBAR/BAR/TI estimators and uncertainty reporting.
- Restraint framework for absolute binding free energy (Boresch-style + standard state corrections).
- Replica scheduling and automatic window equilibration/convergence diagnostics.

### 2) Enhanced sampling
- Metadynamics (well-tempered + multiple walkers).
- Umbrella sampling + WHAM analysis.
- Adaptive biasing force / OPES-style modular bias API.
- Collective variable library (distances, contacts, torsions, RMSD, path CVs).

### 3) Protein–ligand workflow primitives
- Receptor–ligand system builder (protein prep hooks, protonation/tautomer states, cofactors, waters).
- Ligand parameterization adapters (GAFF/OpenFF/CGenFF workflows).
- Pose refinement mini-pipelines (short restrained MD, sidechain relaxation, water equilibration).
- Batch execution for congeneric series (lead optimization mode).

### 4) MM/PBSA and endpoint methods
- Trajectory post-processing for MM/PBSA and MM/GBSA.
- Configurable dielectric models and nonpolar terms.
- Per-residue energy decomposition and interaction hotspot reports.
- Bootstrap confidence intervals for ranking robustness.

### 5) Electrostatics and long-range accuracy
- PME/PPPM implementation for production electrostatics.
- Ewald parameter auto-tuning for target error.
- Long-range LJ corrections and robust tail handling.
- Numerical validation suite for force/energy consistency.

### 6) Quantum-mechanical integrations
- QM/MM interface layer (region definitions, link atoms, electrostatic embedding).
- External QM backend connectors (ORCA, Psi4, xTB) via standardized input/output adapters.
- Semiempirical fast-path for screening and DFT refinement for selected candidates.
- Reaction-coordinate scans and constrained optimizations for mechanism hypotheses.

### 7) Automation, scaling, and reproducibility
- Workflow engine for campaign templates (equilibration, production, analysis, reporting).
- Cloud/HPC schedulers (SLURM/Kubernetes) with restart/checkpoint resilience.
- Artifact tracking (inputs, seeds, environment, hashes) for strict reproducibility.
- Result database + dashboard for portfolio-level decision support.

### 8) AI/ML and Frontier-Hub style integration points
- Feature export for ML reweighting/surrogates (energies, contacts, CV traces).
- Active-learning loop: propose compounds → simulate → retrain prioritization model.
- Uncertainty-aware rank aggregation across physics + ML models.
- API-first service endpoints so design teams can run workflows without MD internals expertise.

### 9) Validation and governance for decision-grade predictions
- Reference benchmark panel (public and internal protein–ligand systems).
- Continuous calibration against experimental ΔG and rank metrics.
- Prediction confidence policies (go/no-go thresholds by project phase).
- Audit-ready reports for model assumptions, parameter provenance, and failure modes.

### Suggested implementation order
1. PME + nonbonded correctness hardening.
2. ABFE/RBFE alchemical core with MBAR/BAR tooling.
3. Enhanced sampling and CV framework.
4. MM/PBSA post-processing.
5. QM/MM adapters.
6. Campaign orchestration, dashboards, and AI integration.

