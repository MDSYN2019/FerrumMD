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
