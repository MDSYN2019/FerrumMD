use crate::lenanrd_jones_simulations::site_site_energy_calculation;
use crate::lennard_jones_simulations::InitOutput; // import the init output struct
use crate::lennard_jones_simulations::Particle; // import the particle struct

use crate::molecule::molecule::System; // import the system struct // input the init struct
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
            // initialize the state with the system configuration and compute the initial potential energy
            let system_potential_energy = 0.0;
            Self {
                system,
                system_potential_energy,
            }
        }

        pub fn update_potential_energy(&mut self) -> () {
            self.system_potential_energy = site_site_energy_calculation(); // TODO
        }
    }
}
