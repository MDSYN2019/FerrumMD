use std::cmp::min;
use crate::lennard_jones_simulations::InitOutput;
use crate::lennard_jones_simulations::Particle; // import the particle struct
use crate::molecule::molecule::System; // import the system struct // input the init struct
use crate::state::state;
use create::constants::BOLTZMANN_CONSTANT;

pub mod replica_exchange {
    pub struct LambdaState {
	pub elec: f64,
	pub vdw: f64,
	pub bonded: f64,
	pub restraints: f64,
    }
    
    pub enum AlchemicalMode {
	Decouple,   // keep intramolecular terms, remove intermolecular
	Annihilate, // remove everything in selected component
    }

    pub struct AlchemicalRegion {
	pub atoms: Vec<usize>,
	pub sterics_mode: AlchemicalMode,
	pub electrostatics_mode: AlchemicalMode,
	pub softcore: SoftcoreParams,
    }

    pub trait Hamiltonian {
	fn energy(&self, state: &SystemState, lambda: &LambdaState) -> f64;
	fn forces(&self, state: &SystemState, lambda: &LambdaState, out: &mut Forces);
	fn du_dlambda(&self, state: &SystemState, lambda: &LambdaState) -> LambdaDerivs;
	fn components(&self, state: &SystemState, lambda: &LambdaState) -> EnergyBreakdown;
    }

    pub struct Replica {
	// represents a single replica in the exchange simulation
	pub id: usize, // unique identifier for the replica
	pub temperature: f64,
	pub state: state,
	pub accepted_exchange: usize,
	pub attempted_exchange: usize,
    }
    
    // methods for the replica, such as updating state, calculating energy
    pub fn attempt_exchange_particles_temperature(replica_1 : &mut Replica, replica_2; &mut Replica) -> () {
	replica_1.attempted_exchange += 1;
	replica_2.attempted_exchange += 2;
	let prob = exchange_probability(replica_1, replica_2);
	if rand::random::<f64>() < prob {
	    // exchange the states of the replicas
	    replica_1.state.swap(&mut replica_2.state);
	    replica_2.state.swap(&mut replica_1.state);	    
	}
    }
    
    pub fn exchange_probability(replica_1 : &mut Replica, replica_2; &mut Replica) -> () {
	beta_1 = 1.0 / (BOLTZMANN_CONSTANT * replica_1.temperature);
	beta_2 = 1.0 / (BOLTZMANN_CONSTANT * replica_2.temperature);
	
	let u_1 = replica1.potential_energy;
	let u_2 = replica2.potential_energy;
	let exponent = (beta_2 - beta_1) * (u_1 - u_2);
	// return the min of 1 and the exponent 
	min(1.0, exponent.exp())
    }
}

