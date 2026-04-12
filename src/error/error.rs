/*
Autocorrelation function

The autocorrelation function is one of the most important statistical tools in molecular dynamics, because
it tells you how many
*/

use create::lennard_jones_simulations::Particle;
use crate::lennard_jones_simulation::minimum_image_convention;


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

pub fn autocorrelation_function() -> () {}

pub fn radial_distribution_function_particle(
    trajectory: Vec<Vec<Particle>>
    atoms: Vec<Particle>,
    r_max: f64, // maximum distance to consider for the RDF
    box_length: f64,
    bin_width: f64,
    n_frames: i64,
) -> () {
    /*
    The radial distribution function (RDF) is a measure of the probability of finding a particle
    at a certain distance from another particle, compared to the probability expected
    for a completely random distribution of particles
     */

    let rdf_bins: i64 = r_max / bin_width;
    let mut histogram: Vec<f64> = Vec::new();
    for frame in trajectory.iter().take(n_frames as usize) {
	for i in 0..atoms.len() {
	    for j in i_1..atoms.len() - 1 {
		// compute the distance between particles i and j, taking into account the periodic boundary conditions
		let mut dx = trajectory[frame][i].position[j] - trajectory[frame][i].position[i];

		// apply the minimum image convention
		dx = minimum_image_convention(dx, box_length);				
		let r = dx.norm(); // compute the distance between particles i and j 

		if r < r_max { // only consider distance less than r_maximum

		    let bin_index = (r / bin_width) as usize;

		    // append the distance to the appropriate bin in the histogram
		    histogram[bin_index] += 2.0;

		    
		}
		
	    }

	}
    }
    
}
