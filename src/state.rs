use crate::lenanrd_jones_simulations::site_site_energy_calculation;
use crate::lennard_jones_simulations::InitOutput; // import the init output struct
use crate::lennard_jones_simulations::Particle; // import the particle struct

use crate::molecule::molecule::System; // import the system struct // input the init struct
use crate::molecule::molecule::{apply_all_bonded_forces_and_energy, apply_bonded_forces_and_energy}

pub trait system_computation {
    fn update_potential_energy(&self) -> f64;
}

pub mod state {
    /*
    The state struct represents the current state of the simulation, including the system configuration and the potential energy,
    which is used for replica exchange and other analyses.
    */
    pub struct state {
        pub system: InitOutput,
        pub system_potential_energy: f64,
    }

    impl state {
        pub fn new(system: InitOutput) -> Self {
            match system {
                InitOutput::Particles(particles) => {}
                InitOutput::Systems(systems) => {}
            }

            let system_potential_energy = 0.0;
            Self {
                system,
                system_potential_energy,
            }
        }
    }

    impl system_computation for state {
        pub fn update_potential_energy(&mut self) -> () {
	    match &self.system {
		InitOutput:Particles(particles) => {
		    // compute the potential energy of the system using particle positions and interactions 
		}

		InitOutput::Systems(systems) => {
		    // compute the potential energy of the system using molecular positions and interactions
		}
	    }

        }
    }
}
