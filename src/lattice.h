#pragma once

/* 

Paper for reference to n-fold algorithm & equations:
https://digitalcommons.kettering.edu/cgi/viewcontent.cgi?article=1010&context=physics_facultypubs

A Fast Serial Algorithm for the Finite Temperature Quenched Potts Model (1993)
Gregory N. Hassold, Elizabeth A. Holm

Paper for reference to probability and energy equations:
http://mimp.materials.cmu.edu/rohrer/theses/151219_final_WilliamFrazierThesis.pdf

A Potts Model Investigation of Complexion Transitions and Abnormal Grain Growth (2015)
William E. Frazier

*/

#include "types.h"
#include "voxel.h"
#include "octree3.h"
#include "boundaries2.h"

#include <cmath>
#include <random>
#include <set>
#include <fstream>
#include <algorithm>

// An object representing a voxel lattice.
class lattice_t
{
private:
	// X, Y, and Z offsets for each neighbor 0-25.
	// NEIGHBOR_LOOKUP tables allow for traversal of voxel neighbors within a one-dimensional loop.
	char NEIGHBOR_LOOKUP_X[NEIGH_COUNT];
	char NEIGHBOR_LOOKUP_Y[NEIGH_COUNT];
	char NEIGHBOR_LOOKUP_Z[NEIGH_COUNT];

	// A lookup table for the e^(-dE / kT) term for each possible dE (-26 to 26) since repeated computation may be expensive.
	activ_t PROB_ETERM_LOOKUP[NEIGH_COUNT * 2 + 1];

	// Temperature that the simulation should run at.
	const activ_t kT = 0.5;

	// Build the neighbor offset and e-term lookup tables.
	void build_lookup_tables()
	{
		std::cout << "Building look up tables..." << std::endl;
		// "Neighbors" of a voxel are all 26 surrounding voxels.
		// Neighbors on a corner have the same "neighborness" as a neighbor sharing a face, even if technically further away.
		// Don't ask me if this is correct... this is how Holm's Fortran code did it.
		char offset_index = 0;
		for (char z = -1; z <= 1; ++z)
			for (char y = -1; y <= 1; ++y)
				for (char x = -1; x <= 1; ++x)
				{
					if (x == 0 && y == 0 && z == 0) continue;

					NEIGHBOR_LOOKUP_X[offset_index] = x;
					NEIGHBOR_LOOKUP_Y[offset_index] = y;
					NEIGHBOR_LOOKUP_Z[offset_index] = z;
					++offset_index;
				}
		
		for (char de = -NEIGH_COUNT; de <= NEIGH_COUNT; ++de)
		{
			PROB_ETERM_LOOKUP[de + NEIGH_COUNT] = exp((-de) / (kT));
		}
		std::cout << "Look up table x: " << NEIGHBOR_LOOKUP_X << std::endl;
		std::cout << "Look up table y: " << NEIGHBOR_LOOKUP_X << std::endl;
		std::cout << "Look up table z: " << NEIGHBOR_LOOKUP_X << std::endl;
		//sleep(5);
	}

	// Get the mobility between two grains.
	activ_t get_mobility(spin_t a, spin_t b)
	{
		return boundary_tracker.is_transformed(a, b) ? transitioned_mobility : default_mobility;
	}

	const char NO_DELTA_E_NEIGHBOR = -50;
	// Calculate the change in energy associated with flipping a voxel to a new spin.
	char get_deltaE(coord_t x, coord_t y, coord_t z, size_t new_spin)
	{
		spin_t curr_spin = voxel_at(x, y, z)->spin,
			   nspin;
		char output = 0;
		bool found = false;

		// dE is equal to the number of neighboring voxels with the current spin minus the number of neighboring voxels with the new spin.
		for (char n = 0; n < NEIGH_COUNT; ++n)
		{
			nspin = neighbor_at(x, y, z, n)->spin;
			if (nspin == new_spin)
			{
				--output;
				found = true;
			}
			else if (nspin == curr_spin)
			{
				++output;
			}
		}

		return found ? output : NO_DELTA_E_NEIGHBOR;
	}

	// Get the probability of a voxel flipping to a new spin.
	// Calculated from Eq. 4.2 on page 42 of Frazier PhD thesis.
	activ_t get_prob(coord_t x, coord_t y, coord_t z, size_t new_spin)
	{
		spin_t curr_spin = voxel_at(x, y, z)->spin;
		if (new_spin == curr_spin) return 0;

		char dE = get_deltaE(x, y, z, new_spin);
		if (dE == NO_DELTA_E_NEIGHBOR) return 0;
		else if (dE < 0) return get_mobility(curr_spin, new_spin);
		else return get_mobility(curr_spin, new_spin) * PROB_ETERM_LOOKUP[dE + NEIGH_COUNT];
	}

	// Get a random float value between min and max.
	activ_t rng(activ_t min, activ_t max)
	{
		//return (((activ_t)rand() / (activ_t)RAND_MAX) * (max - min)) + min;
		return (rng_dis(rng_gen) * (max - min)) + min;
	}

	// Clear and recalculate the overall activity for a voxel.
	void rebuild_voxel_activity(coord_t x, coord_t y, coord_t z)
	{
		//std::cout << "Rebuilding voxel activity..." << std::endl;
		voxel_t *v = voxel_at(x, y, z);
		// Expand on voxel operations in comments.!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		for (char n = 0; n < NEIGH_COUNT; ++n)
		{
			size_t nspin = neighbor_at(x, y, z, n)->spin;
			if (nspin == v->spin || v->has_neighbor(nspin)) continue;

			activ_tree->delta(x, y, z, v->set_neighbor(nspin, get_prob(x, y, z, nspin), &boundary_tracker));
		}
	}
	// Recalculate the activity for a single neighbor of a voxel.
	void rebuild_neighbor_activity(coord_t x, coord_t y, coord_t z, size_t nspin)
	{
		x = (x + side_length_x) % side_length_x;
		y = (y + side_length_y) % side_length_y;
		z = (z + side_length_z) % side_length_z;

		voxel_t *v = voxel_at(x, y, z);

		activ_t new_prob = get_prob(x, y, z, nspin);
		activ_tree->delta(x, y, z, v->set_neighbor(nspin, new_prob, &boundary_tracker));
	}

	// Flip a voxel to a new spin.
	// NOTE: Due to the fact that neighboring spins are accessed/updated, this prevents the simulation from being easily parallelizable (among many other things).
	void flip_voxel(coord_t x, coord_t y, coord_t z, spin_t new_spin)
	{
		voxel_t *v = voxel_at(x, y, z);
		spin_t old_spin = v->spin;
		activ_tree->delta(x, y, z, v->reset(&boundary_tracker));
		v->spin = new_spin;

		rebuild_voxel_activity(x, y, z);
		for (char n = 0; n < NEIGH_COUNT; ++n)
		{
			rebuild_neighbor_activity(x + NEIGHBOR_LOOKUP_X[n], y + NEIGHBOR_LOOKUP_Y[n], z + NEIGHBOR_LOOKUP_Z[n], old_spin);
			rebuild_neighbor_activity(x + NEIGHBOR_LOOKUP_X[n], y + NEIGHBOR_LOOKUP_Y[n], z + NEIGHBOR_LOOKUP_Z[n], new_spin);
		}

		boundary_tracker.track_flip(old_spin, new_spin);

		++total_flips;
		if (boundary_tracker.is_transformed(old_spin, new_spin)) ++transformed_flips;
	}

	// Probabilistically find a voxel that can be flipped based on a random activity (1..system_activity).
	void find_voxel(activ_t desired_activ, coord_t *outx, coord_t *outy, coord_t *outz)
	{
		activ_tree->get_voxel_from_sum_activity(outx, outy, outz, desired_activ, voxels, side_length_x, side_length_y, side_length_z);
	}

	void from_index(size_t index, coord_t *outx, coord_t *outy, coord_t *outz)
	{
		size_t temp;

		temp = index % side_length_x;
		*outx = temp;
		index -= temp;
		index /= side_length_x;

		temp = index % side_length_y;
		*outy = temp;
		index -= temp;
		index /= side_length_y;

		temp = index % side_length_z;
		*outz = temp;
	}

	std::mt19937 rng_gen;
	std::uniform_real_distribution<> rng_dis;

public:
	// The length of x side of the lattice.
	coord_t side_length_x;
	coord_t side_length_y;
	coord_t side_length_z;
	//coord_t side_length;
	int z_propagation_plane;

	// An array of all voxels within the lattice.
	voxel_t *voxels;
	// A counter on the total number of flips the simulation has conducted so far.
	size_t total_flips;
	octree3_t *activ_tree;
	activ_t default_mobility, transitioned_mobility;
	size_t transformed_flips;

	// Total number of possible grains in the simulation.
	// In Holm's code this is a constant value, here we set it to the number of grains within the initial state.
	spin_t grain_count;

	// An object that tracks and controls grain boundary transformations.
	boundary_tracker_t boundary_tracker;

	// Get the overall activity within the lattice.
	activ_t system_activity()
	{
		return activ_tree->system_activity();
	}

	// Constructor for a lattice object.
	lattice_t(coord_t dim_x, coord_t dim_y, coord_t dim_z)
	{
		side_length_x = dim_x;
		side_length_y = dim_y;
		side_length_z = dim_z;

		voxels = new voxel_t[side_length_x * side_length_y * side_length_z];
		total_flips = 0;

		default_mobility = 0.002;
		transitioned_mobility = 0.04;

		total_flips = transformed_flips = 0;

		// the z_plane that should be set for propagation
		z_propagation_plane = 0;

		// Make the area managed by the octree have a side length equal to the next power of two after the real side length.
		// The purpose is to prevent unexpected behavior from integer division (which I spent hours trying to debug...).
		// It results in more memory usage, but (hopefully) shouldn't be a huge problem.
		size_t next_highest_power_of_2 = 1;
		// update: find the largest side length
		coord_t max_side_length = std::max({side_length_x, side_length_y, side_length_z});
		std::cout << "Max side length: " << max_side_length << std::endl;
		while (next_highest_power_of_2 < max_side_length)
		{
			next_highest_power_of_2 *= 2;
		}
		activ_tree = new octree3_t(next_highest_power_of_2, log2(next_highest_power_of_2) + 1);

		rng_gen = std::mt19937(1337);
		rng_dis = std::uniform_real_distribution<>(0.0, 1.0);

		std::cout << "Created lattice of size " << side_length_x << ", " << side_length_y << ", " << side_length_z  << std::endl;
	}
	~lattice_t()
	{
		delete voxels;
		delete activ_tree;
	}

private:
	bool log_transitions;
	std::ofstream transition_log_file;
	double log_timestep;

public:
	// Start logging transitions to the log file.
	void begin_logging_transitions(std::string output_folder)
	{
		std::cout << "Starting to log transitions..." << std::endl;

		output_folder += "transitions.txt";
		log_transitions = true;
		transition_log_file = std::ofstream(output_folder.c_str());
	}
	// Stop logging transitions to the log file.
	void stop_logging_transitions()
	{
		log_transitions = false;
		transition_log_file.close();
	}
	// Write changes to the log file.
	void flush_log_file()
	{
		std::flush(transition_log_file);
	}
	// Set the log's current timestep.
	void set_log_timestep(double timestep)
	{
		log_timestep = timestep;
	}

	// Get the voxel index at the given coordinates (wraps).
	size_t index_at(coord_t x, coord_t y, coord_t z)
	{
		x = (x + side_length_x) % side_length_x;
		y = (y + side_length_y) % side_length_y;
		z = (z + side_length_z) % side_length_z;

		return (x + (y * side_length_x) + (z * side_length_x * side_length_y));
	}
	// Get the voxel at the given coordinates (wraps).
	voxel_t *voxel_at(coord_t x, coord_t y, coord_t z)
	{
		x = (x + side_length_x) % side_length_x;
		y = (y + side_length_y) % side_length_y;
		z = (z + side_length_z) % side_length_z;

		return &voxels[x + (y * side_length_x) + (z * side_length_x * side_length_y)];
	}
	// Get the n'th neighbor of the voxel at the given coordinates.
	voxel_t *neighbor_at(coord_t x, coord_t y, coord_t z, char n)
	{
		return voxel_at(x + NEIGHBOR_LOOKUP_X[n], y + NEIGHBOR_LOOKUP_Y[n], z + NEIGHBOR_LOOKUP_Z[n]);
	}

	// Initialize the lattice (used to build initial activity values at the start of the simulation).
	void init()
	{
		std::cout << "Initializing..." << std::endl;

		build_lookup_tables();

		std::unordered_set<spin_t> spins;

		for(coord_t z = 0; z < side_length_z; ++z)
			for (coord_t y = 0; y < side_length_y; ++y)
				for (coord_t x = 0; x < side_length_x; ++x)
				{
					voxel_at(x, y, z)->index = index_at(x, y, z);
					if (grain_count <= 0) spins.insert(voxel_at(x, y, z)->spin);

					rebuild_voxel_activity(x, y, z);
				}

		if (grain_count <= 0)
		{
			grain_count = spins.size();
		}
		std::cout << "Grain count: " <<  grain_count << std::endl;
		std::cout << "spins count: " <<  spins.size() << std::endl;
		std::cout << "Done initializing." << std::endl;
	}

	// Step the simulation forward, performing a single voxel flip (returns the number of timesteps that the flip theoretically took).
	double step()
	{
		activ_t rand_activ;
		// Sometimes rand_activ is greater than system_activity() (I think), so this is a hack to prevent that.
		do
		{
			rand_activ = rng(0, system_activity());

		} while (rand_activ >= system_activity());

		coord_t vx, vy, vz;
		find_voxel(rand_activ, &vx, &vy, &vz);

		if (!voxel_at(vx, vy, vz)->activity)
		{
			std::cout << "ERROR: Chose a 0-activity voxel. Exiting..." << std::endl;
			exit(0);
		}

		do
		{
			rand_activ = rng(0, voxel_at(vx, vy, vz)->activity);

		} while (rand_activ >= voxel_at(vx, vy, vz)->activity);

		spin_t new_spin = voxel_at(vx, vy, vz)->choose_neighbor(rand_activ);
		flip_voxel(vx, vy, vz, new_spin);

		// This expression is taken from Eq. 20 in Hassold/Holm 1993.
		//return -((double)grain_count - 1) * log(rng(0.01f, 0.99f)) / system_activity();
		// sub for direct compare
		return -log(rng(0.01f, 0.99f)) / system_activity();
	}

private:
	void transition_boundary(boundary_t *boundary)
	{
		//cross_zplane(boundary);
		boundary_tracker.mark_transformed(boundary);
		//std::cout << "Boundaries between " << boundary->a_spin << " and " << boundary->b_spin << std::endl;
		// Update voxel activities and octree for all voxels on the boundary.
		for (auto bvox_iter = boundary->boundary_voxel_indices.begin(); bvox_iter != boundary->boundary_voxel_indices.end(); ++bvox_iter)
		{
			coord_t x, y, z;
			from_index(*bvox_iter, &x, &y, &z);
			rebuild_voxel_activity(x, y, z);
			for (char n = 0; n < NEIGH_COUNT; ++n)
			{
				rebuild_neighbor_activity(x + NEIGHBOR_LOOKUP_X[n], y + NEIGHBOR_LOOKUP_Y[n], z + NEIGHBOR_LOOKUP_Z[n], voxels[*bvox_iter].spin);
			}
		}

		if (log_transitions)
		{
			transition_log_file << boundary->a_spin << '\t' << boundary->b_spin << '\t' << std::to_string(log_timestep) << '\n';
		}
	}

public:
	// Add-on: Check if a boundary cross the z-plane.
	// z_tolerance is for floating points errror
	// It needs: 1. a pair of adjacent voxels from different grain
	// 2. the pair is adjacent to z-boundary
	// revision 1/28/26, adjust to customized z-plane, value should be cfg.z_plane
	bool cross_zplane(boundary_t *boundary ,int Z_TOLERANCE=1)
	{
		//std::cout << "Boundaries between " << boundary->a_spin << " and " << boundary->b_spin << std::endl;
		// Filter z coordinates in the set by z-plane.
		// filter_a holds the voxel that's on the boundary
		std::unordered_set<int> filtered_a, filtered_b;
		for (auto bvox_iter = boundary->boundary_voxel_indices.begin(); bvox_iter != boundary->boundary_voxel_indices.end(); ++bvox_iter)
		{
			coord_t x, y, z;
			from_index(*bvox_iter, &x, &y, &z);
			// Case 1: on both ends, since it wraps around there will be two ends
			if (z_propagation_plane == 0){
				if(abs(z - side_length_z) <= Z_TOLERANCE || abs(z-0) <= Z_TOLERANCE){
					//std::cout << "x, y, z voxel coordinates: " << x  <<", "  << y<<", " << z << std::endl;
					//std::cout << "Current voxel index: " << index_at(x,y,z) << " and spin index : " << voxel_at(x, y, z)->spin <<std::endl;
					if (voxel_at(x, y, z)->spin == boundary->a_spin) {
						filtered_a.insert(*bvox_iter);
					}else{
						filtered_b.insert(*bvox_iter);
					}
				}
			}
			// Case 2: anything in the middle of the grain should not wrap around
			else{
				if(abs(z-z_propagation_plane) <= Z_TOLERANCE){
					if (voxel_at(x, y, z)->spin == boundary->a_spin) {
						filtered_a.insert(*bvox_iter);
					}else{
						filtered_b.insert(*bvox_iter);
					}
				}
			}
			
		}
		// Phase 2:
		// for each voxel in the blue set, for each of its neighbor off set
		for (auto filtered_vox : filtered_a) 
		{
			coord_t fx, fy, fz;
			from_index(filtered_vox, &fx, &fy, &fz);
			for (int z = -1; z <= 1; ++z)
				for (int y = -1; y <= 1; ++y)
					for (int x = -1; x <= 1; ++x)
					{
						if (x == 0 && y == 0 && z == 0) continue;
						if (filtered_b.count(index_at(fx+x, fy+y, fz+z))){
							std::cout << "[INFO] Boundaries between " << boundary->a_spin << " and " << boundary->b_spin << std::endl;
							std::cout << "[INFO] x, y, z voxel coordinates: " <<fx  <<", "  << fy<<", " << fz << std::endl; 
							std::cout << "[INFO] Found: " << fx+x  <<", "<<  fy+y <<", " << fz+z << std::endl;
							return true;
						}
					}
		}
		return false;

		// Filter z coordinates in the set by z-plane.
		// get unordered_set for each two grains. call red and blue
		// for each voxel in the blue set, for each of its neighbor off set
		// find if the red set contains the neighbor. If found, return true.
	}

	// Transitition a certain number of random grain boundaries.
	// Yes, I know this function is a mess. It can probably be simplified a bit...
	void transition_boundaries(size_t count, double propagation_chance, double propagation_ratio, bool use_potential_energy)
	{
		std::cout << "Transitioning " << count << " boundaries..." << std::endl;

		std::set<size_t> flip_indices;
		std::set<size_t> propagate_indices;
		// what does this do?????
		boundary_tracker.remove_bad_boundaries();

		if (count > boundary_tracker.total_boundary_count - boundary_tracker.transformed_boundary_count)
		{
			count = boundary_tracker.total_boundary_count - boundary_tracker.transformed_boundary_count;
		}

		size_t
			propagate_count = count * propagation_chance,
			flip_count = count - propagate_count;

		if (boundary_tracker.transformed_boundary_count < propagate_count)
		{
			propagate_count = boundary_tracker.transformed_boundary_count;
			flip_count = count - propagate_count;
		}
		
		// flip_indices and propagate_indices are ordered sets that contain the "indices" of untransformed and
		// transformed boundaries to flip/propagate on. Since unordered maps cannot be easily indexed, and
		// the overall collections of untransformed and transformed boundaries are combined, each "index" refers
		// to the i'th instance of that type of boundary. For example, if flip_indices contains the value "3", that
		// tells the algorithm to flip the 4th (because of zero-indexing) untransformed boundary that it comes across.
		size_t index;
		for (size_t i = 0; i < flip_count; ++i)
		{
			do
			{
				index = rng(0, boundary_tracker.total_boundary_count - boundary_tracker.transformed_boundary_count);
			} while (flip_indices.find(index) != flip_indices.end());
			flip_indices.insert(index);
		}
		// debug purpose
		// std::cout << "Flip indices: " << std::endl;
		// for (auto flip_index : flip_indices) {
		// 	std::cout << "index: " << flip_index  << "   ";
		// 	// convert the index to x,y, and z coordinates. Check if the z boundary modify the side length.
		// 	// check if the flip pass the z boundary. Convert the flip index to boundary

		// }
		// std::cout << "\n Total flip indices: " << flip_indices.size() << std::endl;

		for (size_t i = 0; i < propagate_count; ++i)
		{
			do
			{
				index = rng(0, boundary_tracker.transformed_boundary_count);
			} while (propagate_indices.find(index) != propagate_indices.end());
			propagate_indices.insert(index);
		}
		// debug purpose
		// std::cout << "Propagate indices: " << std::endl;
		// for (auto propagate_index : propagate_indices) {
		// 	std::cout << "index: " << propagate_index  << "   "; 
		// 	coord_t x, y, z;
		// 	from_index(propagate_index, &x, &y, &z);
		// 	std::cout << "x: " << x << ", y: " <<y  << ", z: "<< z  << "   \n"; 
		// 	// convert the index to x,y, and z coordinates. Check if the z boundary modify the side length.
		// }
		// std::cout << "\n Total propagate indices: " << propagate_indices.size() << std::endl;
		size_t num_random_propagated = 0, num_random_flipped = 0, num_poteng_propagated = 0;

		bool end_loop = false;
		std::set<size_t>::iterator flip_iter = flip_indices.begin(), propagate_iter = propagate_indices.begin();

		// These count the number of transformed and untransformed boundaries come across so far.
		size_t untrans_count = 0, trans_count = 0;

		// Iterate over the entire boundary map.
		for (auto sm_iter = boundary_tracker.boundary_map.begin(); sm_iter != boundary_tracker.boundary_map.end(); ++sm_iter)
		{
			// debug
			//std::cout << "Sm_iter: " << &sm_iter << std::endl;
			//std::cout << "Untrans_count: " << untrans_count << std::endl;
			for (auto lg_iter = sm_iter->second.begin(); lg_iter != sm_iter->second.end(); ++lg_iter)
			{
				//std::cout << "lg_iter: " << lg_iter << std::endl;
				boundary_t *boundary = lg_iter->second;
				// If the boundary is transformed, check if we should propagate from it.
				if (boundary->transformed)
				{
					if (propagate_iter != propagate_indices.end())
					{
						if (*propagate_iter == trans_count)
						{
							size_t prop_num = boundary->junctions.size() * propagation_ratio;
							if (propagation_ratio <= 0) prop_num = 1;

							// !!! This does not choose a random junction, but just goes in order.
							bool found_junc = false;
							for (auto junc_iter = boundary->junctions.begin(); junc_iter != boundary->junctions.end(); ++junc_iter)
							{
								if (!junc_iter->first->transformed)
								{
									transition_boundary(junc_iter->first);

									found_junc = true;
									++num_random_propagated;

									// This is a hack to early-stop propagation. It should be cleaner but I am very tired...
									if (num_random_propagated >= propagate_count)
									{
										propagate_iter = propagate_indices.end();
										--propagate_iter;
										break;
									}

									--prop_num;
									if(prop_num <= 0) break;
								}
							}

							if (!found_junc)
							{
								// This could cause problems if trans_count = total transformed boundary count...
								size_t new_prop_index = trans_count;
								while (propagate_indices.find(++new_prop_index) != propagate_indices.end());
								propagate_indices.insert(new_prop_index);
							}

							++propagate_iter;
						}
						++trans_count;
					}

					if (use_potential_energy)
					{
						// Potential energy propagation.
						if (boundary->previous_surface_area != 0)
						{
							boundary->potential_energy += boundary->previous_surface_area - boundary->area();
							if (boundary->potential_energy < 0)
							{
								boundary->potential_energy = 0;
							}

							bool potential_propagation = true;
							while (potential_propagation)
							{
								potential_propagation = false;
								boundary_t *smallest_junc = nullptr;

								for (auto junc_iter = boundary->junctions.begin(); junc_iter != boundary->junctions.end(); ++junc_iter)
								{
									if (!junc_iter->first->transformed)
									{
										if (smallest_junc == nullptr || smallest_junc->area() > junc_iter->first->area()) smallest_junc = junc_iter->first;
									}
								}

								if (smallest_junc != nullptr && smallest_junc->area() <= boundary->potential_energy)
								{
									transition_boundary(smallest_junc);
									boundary->potential_energy -= smallest_junc->area();
									potential_propagation = true;
									++num_poteng_propagated;
								}
							}
						}
						boundary->previous_surface_area = boundary->area();
					}
				}
				// If the boundary is untransformed, check if we should flip it.
				else if (!boundary->transformed && flip_iter != flip_indices.end())
				{
					// REVISION 1/28, test for different planes
					if (!cross_zplane(boundary)) {
						continue;
					}
					if (*flip_iter == untrans_count)
					{
						//std::cout << "Flip iter: " << *flip_iter << std::endl;
						transition_boundary(boundary);

						++num_random_flipped;
						++flip_iter;
					}
					++untrans_count;
				}
				// If we have exhausted the flip and propagate lists, end the loop.
				else if (flip_iter == flip_indices.end() && propagate_iter == propagate_indices.end())
				{
					end_loop = true;
					break;
				}
			}
			if (end_loop) break;
		}

		std::cout << "Transitioned boundaries: " << boundary_tracker.transformed_boundary_count << " / " << boundary_tracker.total_boundary_count << " boundaries..." << std::endl;
		std::cout << "# Transitioned via propagation: " << num_random_propagated << ", via random flipping: " << num_random_flipped << ", via potential energy: " << num_poteng_propagated << "..." << std::endl;
	}
};