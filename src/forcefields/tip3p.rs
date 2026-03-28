use nalgebra::Vector3;
// Oxygen and Hydrogen charge
const Q_O: f64 = -0.834;
const Q_H: f64 = 0.417;

// LJ parameters
const SIGMA: f64 = 3.1507;
const EPSILON: f64 = 0.1521;

struct tip3p {
    oxygen: Vector3<f64>,
    hydrogen1: Vector3<f64>,
    hydrogen2: Vector3<f64>,
}
