
use crate::lennard_jones_simulations::Particle;
use crate::molecule::molecule::System;
use nalgebra::Vector3;
/*
Create the subcells required for efficiently computing the intermolecular interactions
within a certain radius cutoff rather than working with all the atoms in the system
 */

fn wrap_index(i: isize, n: usize) -> usize {
    let n = n as isize;
    (((i % n) + n) % n) as usize
}

fn cell_id(ix: usize, iy: usize, iz: usize, n: usize) -> usize {
    (ix * n + iy) * n + iz
}

//fn cell_id_to_3d

fn position_to_cell_3d(
    pos: &Vector3<f64>,
    box_size: &Vector3<f64>,
    n_cells: usize,
) -> (usize, usize, usize) {
    /*
    Map to [0, L) first if you're storing unwrapped coordinates
    if your pos is already in [0, L), you can skip this
     */

    let x = pos.x.rem_euclid(box_size.x); // wrap around x coordinate
    let y = pos.y.rem_euclid(box_size.y); // wrap around y coordinate
    let z = pos.z.rem_euclid(box_size.z); // wrap around z coordinate

    // normalized coordinates around [0, 1)

    let fx = (x / box_size.x) * n_cells as f64;
    let fy = (y / box_size.y) * n_cells as f64;
    let fz = (z / box_size.z) * n_cells as f64;

    // floor gives 0..n_cells-1, but we still wrap to be safe
    let ix = wrap_index(fx.floor() as isize, n_cells);
    let iy = wrap_index(fy.floor() as isize, n_cells);
    let iz = wrap_index(fz.floor() as isize, n_cells);

    (ix, iy, iz)
}

pub struct SimulationBox {
    pub x_dimension: f64,
    pub y_dimension: f64,
    pub z_dimension: f64,
}

pub struct MolecularCoordinates {
    /*
    struct to store information about each subdivsion of cells
    */
    pub center: Vector3<f64>,
    //pub length: f64,
    pub half_length: Vector3<f64>,
    pub index: Vector3<usize>,
    pub atom_index: Vec<usize>,
    pub head: Vec<Option<usize>>,
    pub next: Vec<Option<usize>>,
}

impl SimulationBox {
    pub fn create_subcells(&self, n_cells: usize) -> Vec<MolecularCoordinates> {
        let n_cells_per_dim = n_cells;
        let n_cells_total = n_cells_per_dim * n_cells_per_dim * n_cells_per_dim;
        let mut cells = Vec::with_capacity(n_cells_total); // create the cells

        let dx = self.x_dimension / n_cells_per_dim as f64;
        let dy = self.y_dimension / n_cells_per_dim as f64;
        let dz = self.z_dimension / n_cells_per_dim as f64;

        for ix in 0..n_cells_per_dim {
            for iy in 0..n_cells_per_dim {
                for iz in 0..n_cells_per_dim {
                    // create the cells and push the coordinates (with the center implemented)
                    cells.push(MolecularCoordinates {
                        center: Vector3::new(
                            (ix as f64 + 0.5) * dx,
                            (iy as f64 + 0.5) * dy,
                            (iz as f64 + 0.5) * dz,
                        ),
                        half_length: Vector3::new(dx * 0.5, dy * 0.5, dz * 0.5),
                        index: Vector3::new(ix, iy, iz),
                        atom_index: Vec::new(),
                        head: vec![None; n_cells_total],
                        next: Vec::new(), //
                    });
                }
            }
        }

        cells
    }

    pub fn store_atoms_in_cells_particles(
        &self,
        particles: &mut Vec<Particle>,
        cells: &mut Vec<MolecularCoordinates>, // the created cells
        n_cells: usize,
    ) -> () {
        /*
        I don't need to store all the coordinates, just the indices of the atoms/systems
         */

        // Clear previous contents
        for cell in cells.iter_mut() {
            cell.atom_index.clear();
        }

        let box_size = Vector3::new(self.x_dimension, self.y_dimension, self.z_dimension);
        // let n_cells_total = cells.len();

        // Store all the particles in the appropriate

        for (i, particle) in particles.iter_mut().enumerate() {
            let (ix, iy, iz) = position_to_cell_3d(&particle.position, &box_size, n_cells);

            let cid = cell_id(ix, iy, iz, n_cells);
            cells[cid].atom_index.push(i);
        }
    }

    pub fn store_atoms_in_cells_systems(
        &self,
        systems: &mut Vec<System>,
        cells: &mut Vec<MolecularCoordinates>,
        n_cells: usize,
    ) -> () {
        /*
        I don't need to store all the coordinates, just the indices of the atoms/systems
         */

        // Clear previous contents
        for cell in cells.iter_mut() {
            cell.atom_index.clear();
        }

        let box_size = Vector3::new(self.x_dimension, self.y_dimension, self.z_dimension);

        for (i, system) in systems.iter_mut().enumerate() {
            for particle in system.atoms.iter() {
                let (ix, iy, iz) = position_to_cell_3d(&particle.position, &box_size, n_cells);
                let cid = cell_id(ix, iy, iz, n_cells);
                cells[cid].atom_index.push(i);
            }
        }
    }
}
