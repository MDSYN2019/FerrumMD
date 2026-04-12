/*
Autocorrelation function

The autocorrelation function is one of the most important statistical tools in molecular dynamics, because
it tells you how many
*/

use crate::lennard_jones_simulations::minimum_image_convention;
use crate::lennard_jones_simulations::Particle;

struct rdf_results {
    hist: Vec<f64>,
    g: Vec<f64>,
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
    trajectory: Vec<Vec<Particle>>,
    atoms: Vec<Particle>,
    r_max: f64, // maximum distance to consider for the RDF
    box_length: f64,
    bin_width: f64,
    n_frames: i64,
) -> rdf_results {
    /*
    The radial distribution function (RDF) is a measure of the probability of finding a particle
    at a certain distance from another particle, compared to the probability expected
    for a completely random distribution of particles

    The radial distribution function (RDF) is a measure of the probability of findind a particle
    at a certain distance from another particle, compared to the probability expected for a completely random
    distribution of particles
     */

    let rdf_bins: f64 = r_max / bin_width;
    let mut histogram = vec![0.0; rdf_bins as usize]; // initialize the histogram with the number of bins
    let mut g = vec![0.0; rdf_bins as usize]; // initialize the RDF historgram with the number of bings

    for frame in trajectory.iter() {
        for i in 0..atoms.len() {
            for j in i..atoms.len() - 1 {
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
    let mut volume: f64 = box_length.powi(3);
    let mut density: f64 = atoms.len() as f64 / volume;

    for bin_index in 0..rdf_bins as usize {
        let mut r_inner: f64 = bin_index as f64 * bin_width; // the inner radius of the bin
        let mut r_outer: f64 = (bin_index as f64 + 1.0) * bin_width;
        let mut r_mid: f64 = r_inner + 0.5 * bin_width; // the mid point of the bin

        let shell_volume: f64 =
            (4.0 / 3.0) * std::f64::consts::PI * (r_outer.powi(3) - r_inner.powi(3));
        let ideal_count = density * shell_volume * atoms.len() as f64;

        g[bin_index] = histogram[bin_index] / ideal_count;
    }

    let result = rdf_results {
        hist: histogram,
        g: g,
    };

    result
}

pub fn autocorrelation_function() -> () {}
