/*
Autocorrelation function

The autocorrelation function is one of the most important statistical tools in molecular dynamics, because
it tells you how many
*/

use crate::lennard_jones_simulations::minimum_image_convention;
use crate::lennard_jones_simulations::Particle;
use std::collections::HashMap;

pub struct RdfResults {
    pub hist: Vec<f64>,
    pub g: Vec<f64>,
}

pub fn compute_average_val(
    container_value: &mut Vec<f32>,
    block_steps: u64,
    number_of_steps: u64,
) -> ()
/*
    Declare the val that we wish to block average; we need to first store the values
    in a dynamic storage (here a vec!) and we then average it through the vec value
    
 */
{
    if block_steps == 0 {
        log::warn!("Invalid block configuration: block_steps must be > 0");
        return;
    }
    if container_value.is_empty() {
        log::warn!("No values provided for block averaging");
        return;
    }
    if number_of_steps % block_steps != 0 {
        log::warn!("Invalid block configuration: number_of_steps={number_of_steps}, block_steps={block_steps}");
    }
    // ensure that the property we want to compute the block average for
    // actually exists within the simulation code

    for (i, chunk) in container_value.chunks(block_steps as usize).enumerate() {
        let summed_values: f32 = chunk.iter().sum();
        let chunk_len = chunk.len() as f32;
        log::info!(
            "Block {i:>4} | avg={:.6} (block size={block_steps})",
            summed_values / chunk_len
        );
    }
}

pub fn radial_distribution_function_particle(
    trajectory: &[Vec<Particle>],
    r_max: f64, // maximum distance to consider for the RDF
    box_length: f64,
    bin_width: f64,
) -> RdfResults {
    /*
    The radial distribution function (RDF) is a measure of the probability of finding a particle
    at a certain distance from another particle, compared to the probability expected
    for a completely random distribution of particles

    The radial distribution function (RDF) is a measure of the probability of findind a particle
    at a certain distance from another particle, compared to the probability expected for a completely random
    distribution of particles
     */

    if trajectory.is_empty() || r_max <= 0.0 || box_length <= 0.0 || bin_width <= 0.0 {
        return RdfResults {
            hist: Vec::new(),
            g: Vec::new(),
        };
    }

    let n_atoms = trajectory[0].len();
    if n_atoms < 2 {
        return RdfResults {
            hist: Vec::new(),
            g: Vec::new(),
        };
    }

    if trajectory.iter().any(|frame| frame.len() != n_atoms) {
        log::warn!("Trajectory contains frames with inconsistent atom counts; returning empty RDF");
        return RdfResults {
            hist: Vec::new(),
            g: Vec::new(),
        };
    }

    let n_frames = trajectory.len() as f64;
    let rdf_bins: f64 = r_max / bin_width;
    let mut histogram = vec![0.0; rdf_bins as usize]; // initialize the histogram with the number of bins
    let mut g = vec![0.0; rdf_bins as usize]; // initialize the RDF historgram with the number of bings

    for frame in trajectory.iter() {
        for i in 0..n_atoms {
            for j in (i + 1)..n_atoms {
                // compute the distance between particles i and j, taking into account the periodic boundary conditions
                let mut dx = frame[j].position - frame[i].position;
                // apply the minimum image convention
                dx = minimum_image_convention(dx, box_length);
                let r = dx.norm(); // compute the distance between particles i and j
                if r < r_max {
                    // only consider distance less than r_maximum
                    let bin_index = (r / bin_width) as usize;
                    // append the distance to the appropriate bin in the histogram
                    histogram[bin_index] += 2.0;
                }
            }
        }
    }
    // normalization
    let volume: f64 = box_length.powi(3);
    let density: f64 = n_atoms as f64 / volume;

    for bin_index in 0..rdf_bins as usize {
        let r_inner: f64 = bin_index as f64 * bin_width; // the inner radius of the bin
        let r_outer: f64 = (bin_index as f64 + 1.0) * bin_width;

        let shell_volume: f64 =
            (4.0 / 3.0) * std::f64::consts::PI * (r_outer.powi(3) - r_inner.powi(3));
        let ideal_count = density * shell_volume * n_atoms as f64 * n_frames;

        g[bin_index] = if ideal_count > 0.0 {
            histogram[bin_index] / ideal_count
        } else {
            0.0
        };
    }

    let result = RdfResults {
        hist: histogram,
        g: g,
    };

    result
}

pub fn mean_squared_displacement_particle(
    trajectory: &[Vec<Particle>],
    number_of_frames: i32,
    box_length: f64,
) -> Result<HashMap<i32, f64>, String> {
    /*
    The mean squared displacement (MSD) is a measure of the average distance that particles have moved
    from their initial positions over time. It is defined as the average of the squared displacements
     */

    let mut msd: HashMap<i32, f64> = HashMap::new(); // initialize a hashmap to store the MSD values for different time lags (tau)
    let number_of_particles: usize = trajectory[0].len(); // get the number of particles from the first frame of the trajectory
    for tau in 0..number_of_frames - 1 {
        // iterate over the time lags (tau) from 0 to the total number of frames minus 1
        //
        let mut sum_squared_displacement: f64 = 0.0;
        let mut count: i32 = 0;

        for t0 in 0..number_of_frames - tau {
            // iterate over the starting time t0 over the time lag tau
            for i in 0..number_of_particles {
                /*
                We wish to store the result for different time lags (tau) in order to compute
                the diffusion coefficient from the slop of the MSD curve at long times
                to get the difference in the position from time 0 to time t
                 */
                let mut diff = trajectory[t0 as usize + tau as usize][i].position
                    - trajectory[t0 as usize][i].position;
                // apply the minimum image convention to account for the periodic boundary condiitons
                diff = minimum_image_convention(diff, box_length);
                let square_displacement = diff.norm_squared(); // compute the squared displacement for particle i at time lag tau
                sum_squared_displacement += square_displacement; // accumulate the squared displacement for all the particles and all the starting times t0
                count += 1;
            }
        }
        // compute the average squared displacement for the time lag tau by dividing the sum of squared displacement by the total number of displacemnets (count)
        if tau % 10 == 0 {
            log::info!("Computes MSD for time lag tau = {tau}");
        }
        if count > 0 {
            msd.insert(tau, sum_squared_displacement / count as f64);
        } else {
            msd.insert(tau, 0.0);
        }
        if msd.is_empty() {
            return Err(
                "Error: no valid MSD values computed; check trajectory data and parameters"
                    .to_string(),
            );
        }
    }
    return Ok(msd);
}

pub fn autocorrelation_function() -> () {}
