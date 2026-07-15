//! FerrumMD demo executable.
//!
//! Demonstrates:
//! - reading CHARMM and GROMACS/ITP force-field files
//! - running a small Lennard-Jones MD example
//! - running a simple molecular H2 system

use log::{error, info, warn};
use std::env;

#[cfg(feature = "mpi")]
use mpi::traits::*;

use sang_md::lennard_jones_simulations;
use sang_md::molecule::{charmm, martini, molecule};

const DEFAULT_INTEGRATOR: &str = "velocity_verlet";
const DEFAULT_TEMPERATURE: f64 = 300.0;
const DEFAULT_BOX_LENGTH: f64 = 10.0;
const DEFAULT_TIME_STEP: f64 = 0.0005;
const DEFAULT_CUTOFF: f64 = 10.0;

fn main() {
    init_logging();
    run_forcefield_examples();

    let integrator = parse_integrator_arg();
    let md_mode = integrator_to_md_mode(&integrator);

    #[cfg(feature = "mpi")]
    let universe = mpi::initialize().expect("MPI initialization failed");

    #[cfg(feature = "mpi")]
    let world = universe.world();

    #[cfg(feature = "mpi")]
    {
        if world.rank() == 0 {
            info!("Running MPI-enabled MD example with integrator={integrator}");
        }
        run_lennard_jones_example_mpi(md_mode, &world);
        run_h2_system_example_mpi(md_mode, &world);
    }

    // serial execution if MPI is not enabled
    #[cfg(not(feature = "mpi"))]
    {
        info!("Running single-process MD example with integrator={integrator}");

        run_lennard_jones_example(md_mode);
        run_h2_system_example(md_mode);
    }
}

fn init_logging() {
    // Set up logging with a default level of INFO, configurable via the RUST_LOG environment variable
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
}

fn parse_integrator_arg() -> String {
    env::args()
        .find_map(|arg| arg.strip_prefix("--integrator=").map(str::to_owned))
        .unwrap_or_else(|| DEFAULT_INTEGRATOR.to_string()) // default to velocity_verlet if no argument is provided
}

fn integrator_to_md_mode(integrator: &str) -> &'static str {
    match integrator {
        "velocity_verlet" => "berendsen",
        "monte_carlo" => "monte_carlo",
        other => {
            warn!(
                "Unknown integrator '{other}', falling back to velocity_verlet \
                 with Berendsen thermostat."
            );
            "berendsen"
        }
    }
}

fn run_forcefield_examples() {
    print_charmm_tip3p_example();
    print_itp_tip3p_example();
}

fn print_charmm_tip3p_example() {
    let path = "/home/sang/Desktop/Sandbox/rust/FerrumMD/ff/tip3p_example/tip3p_sample.rtf";

    let forcefield = match charmm::CharmmForceField::read_file_new(path) {
        Ok(forcefield) => forcefield,
        Err(err) => {
            error!("Failed to read CHARMM TIP3P force field from '{path}': {err}");
            return;
        }
    };

    println!(
        "\nCHARMM TIP3P force field
         \nAtom types: {:?}
         \nResidue names: {:?}
         \nBonds: {:?}\n",
        forcefield.atom_types, forcefield.residue_name, forcefield.bonds,
    );
}

fn print_itp_tip3p_example() {
    let path = "ff/Charmm27.ff/charmm27.ff/tip3p.itp";

    let forcefield = match martini::MartiniForceField::read_itp(path) {
        Ok(forcefield) => forcefield,
        Err(err) => {
            error!("Failed to read ITP TIP3P force field from '{path}': {err}");
            return;
        }
    };

    println!(
        "\nITP TIP3P force field
         \nAtoms: {:?}
         \nBonds: {:?}\n",
        forcefield.atoms, forcefield.bonds,
    );
}

fn create_charged_lennard_jones_system() -> Option<sang_md::lennard_jones_simulations::InitOutput> {
    let mut system = match lennard_jones_simulations::create_atoms_with_set_positions_and_velocities(
        3,
        DEFAULT_TEMPERATURE,
        30.0,
        DEFAULT_BOX_LENGTH,
        DEFAULT_BOX_LENGTH,
        false,
    ) {
        Ok(system) => system,
        Err(err) => {
            error!("Failed to create Lennard-Jones atoms: {err}");
            return None;
        }
    };

    if let lennard_jones_simulations::InitOutput::Particles(particles) = &mut system {
        for (idx, particle) in particles.iter_mut().enumerate() {
            particle.charge = if idx % 2 == 0 { 1.0 } else { -1.0 };
        }
    }

    Some(system)
}

#[cfg(not(feature = "mpi"))]
fn run_lennard_jones_example(md_mode: &str) {
    let Some(mut system) = create_charged_lennard_jones_system() else {
        return;
    };

    lennard_jones_simulations::run_md_nve(
        &mut system,
        30,
        DEFAULT_TIME_STEP,
        DEFAULT_CUTOFF,
        md_mode,
        30.0,
    );

    if md_mode != "monte_carlo" {
        lennard_jones_simulations::run_md_nve(
            &mut system,
            3000,
            DEFAULT_TIME_STEP,
            DEFAULT_CUTOFF,
            "andersen",
            30.0,
        );
    }
}

#[cfg(feature = "mpi")]
fn run_lennard_jones_example_mpi(md_mode: &str, world: &impl mpi::traits::Communicator) {
    let Some(mut system) = create_charged_lennard_jones_system() else {
        return;
    };

    lennard_jones_simulations::run_md_nve_mpi(
        &mut system,
        30,
        DEFAULT_TIME_STEP,
        DEFAULT_CUTOFF,
        md_mode,
        world,
    );
}

fn create_h2_systems() -> sang_md::lennard_jones_simulations::InitOutput {
    let h2 = molecule::make_h2_system();
    let mut systems = molecule::create_systems(&h2, 12);

    lennard_jones_simulations::set_molecular_positions_and_velocities(
        &mut systems,
        DEFAULT_TEMPERATURE,
        DEFAULT_BOX_LENGTH,
    );

    systems
}

#[cfg(not(feature = "mpi"))]
fn run_h2_system_example(md_mode: &str) {
    let mut systems = create_h2_systems();

    lennard_jones_simulations::run_md_nve(
        &mut systems,
        30,
        DEFAULT_TIME_STEP,
        DEFAULT_CUTOFF,
        md_mode,
        30.0,
    );
}

#[cfg(feature = "mpi")]
fn run_h2_system_example_mpi(md_mode: &str, world: &impl mpi::traits::Communicator) {
    let mut systems = create_h2_systems();

    lennard_jones_simulations::run_md_nve_mpi(
        &mut systems,
        30,
        DEFAULT_TIME_STEP,
        DEFAULT_CUTOFF,
        md_mode,
        world,
    );
}
