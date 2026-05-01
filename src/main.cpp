// (c) 2023 Benjamin Zalatan Productions
// Extended and modified by Surui Huang, 2025
// Keep circulating the tapes

#include <iostream>
#include <vector>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <zip.h>

#include "vtk.h"
#include "lattice.h"
#include "octree3.h"
#include "debug_timer.h"
#include "config.h"
#include "analysis.h"
#include "zip_utils.h"


// cd C:\Stuff\School\summer 2023\grainsim
// g++ -O3 CPPGrainSim/main.cpp -o grainsim.out -static
void print_usage()
{
    std::cerr << "Usage:\n"
			<< "  grainsim [--config <file>] [--initial <seed.ph>] [--output <folder>]\n\n"
			<< "Options:\n"
			<< "  --config <file>    Path to config file (default: grainsim_config.txt)\n"
			<< "  --initial <path>   Override INITIAL_STATE_FILE from config\n"
			<< "  --output <path>    Override OUTPUT_FOLDER from config\n"
			<< std::endl;
}


int main(int argc, char *argv[])
{
	// add on for params config 
	std::string config_path = "grainsim_config.txt";
	std::string init_override;
	std::string out_override;
	int transition_override = -1;   // -1 = don't override

	for(int i = 1; i < argc; ++i) {
		std::string arg = argv[i];
		// config file path
		if (arg == "--config" && i+1 < argc) {
			config_path = argv[++i];
		}
		else if (arg == "--initial" && i+1 < argc) {
			init_override = argv[++i];
		}
		else if (arg == "--output" && i+1 < argc) {
			out_override = argv[++i];
		}
		// override with transition count
		else if (arg == "--transition-count" && i+1 < argc) {
			transition_override = std::atoi(argv[++i]);
		}
		else {
			std::cerr << "Unknown or incomplete argument: " << arg << "\n\n";
			print_usage();
			std::exit(1);
		}
	}
	// Load the config file.
	config_t cfg;
	std::cout << "Loading: " << cfg.initial_state_path << std::endl;
	cfg.load_config(config_path, init_override, out_override, transition_override);

	std::cout << "ramp_up_enabled=" << cfg.ramp_up_enabled
          << " transition_count=" << cfg.transition_count
          << " transition_time_count_increase=" << cfg.transition_time_count_increase
          << " high_nucleation_rate_limit=" << cfg.high_nucleation_rate_limit
          << " transition_interval=" << cfg.transition_interval
          << "\n";

	// Create the lattice from file.
	lattice_t *cube;

	// If the scale multiplier is not 1, scale the lattice.
	if (cfg.scale_multiplier != 1)
	{
		lattice_t *temp = vtk::from_file(cfg.initial_state_path.c_str(), false);
		cube = vtk::scale_lattice(temp, cfg.scale_multiplier, false);
		delete temp;
	}
	else
	{
		cube = vtk::from_file(cfg.initial_state_path.c_str(), false);
	}

	lattice_analyzer_t analyze;

	cube->default_mobility = cfg.default_mobility;
	cube->transitioned_mobility = cfg.transitioned_mobility;
	cube->grain_count = cfg.const_grain_count;
	cube->z_propagation_plane = cfg.z_propagation_plane;
	cube->init();

	// Generate the checkpoint list.
	std::vector<double> checkpoints;
	if (!cfg.checkpoints.empty())
	{
		cfg.checkpoints_to_vector(&checkpoints);
	}
	unsigned char curr_checkpoint = 0;

	// Start global timer.
	debug_timer_t timer;
	timer.start();

	double timestep = 0, curr_step, log_duration = 0, transition_duration = 0, next_checkpoint = cfg.checkpoint_interval;
	int vtkcount = 0;

	if (cfg.log_transitions) cube->begin_logging_transitions(cfg.output_folder);

	// initialize transition_interval as the config file
	double ramp_up_transition_count  = 0;
	// for rampup updates
	double next_transition_time = cfg.transition_interval;

	// create small zips for every batch size
	std::vector<std::string> small_zips;

	// Main simulation loop.
	while (true)
	{
		// Flip a voxel and store elapsed timesteps.
		curr_step = cube->step();

		static double last = -1.0;
		if (ramp_up_transition_count != last) {
			//std::cout << "t=" << timestep << " ramp=" << ramp_up_transition_count << "\n";
			last = ramp_up_transition_count;
		}

		// Update current timestep.
		timestep += curr_step;
		log_duration += curr_step;
		transition_duration += curr_step;

		// updates for ramp-up nucleaton rate
		if (cfg.ramp_up_enabled) {
			// NucleationRate(t) = Constant * MCS, capped at limit, discrete
			if (timestep >= next_transition_time) {
				ramp_up_transition_count = std::min(
					cfg.high_nucleation_rate_limit, 
					ramp_up_transition_count + cfg.transition_time_count_increase
				);
				std::cout << "Ramp up nucleation rate increased to: " << ramp_up_transition_count << std::endl;

				next_transition_time += cfg.transition_interval;
			}
		} else {
			ramp_up_transition_count = cfg.transition_count;
		}

		// Debug logging.
		if (log_duration >= 20000)
		{
			std::cout << "T = " << timestep << ", dT = " << curr_step << ", A = " << cube->system_activity() << ", Flips = " << cube->total_flips << ", tFlips = " << cube->transformed_flips << ", dTime = " << timer.lap() << " sec, tTime = " << timer.total() << " sec" << std::endl;
			log_duration = 0;
		}

		// Transition some boundaries if applicable.
		// if (transition_duration >= cfg.transition_interval && cfg.transition_count > 0)
		if (transition_duration >= cfg.transition_interval && ramp_up_transition_count > 0)
		{
			std::cout << "Transition some boundaries if applicable." << std::endl;
			if (cfg.log_transitions) cube->set_log_timestep(timestep);

			//cube->transition_boundaries(cfg.transition_count, cfg.propagation_chance, cfg.propagation_ratio, cfg.use_potential_energy);
			cube->transition_boundaries(ramp_up_transition_count, cfg.propagation_chance, cfg.propagation_ratio, cfg.use_potential_energy);
			transition_duration = 0;
		}

		//std::cout << "Current Timestep: " << timestep << std::endl;
		// Check if VTK should be generated.
		if (checkpoints.size() > 0 && curr_checkpoint < checkpoints.size() && timestep >= checkpoints[curr_checkpoint]) // The current timestep is an explicit checkpoint.
		{
			// WRITE TO VTK 
			std::cout << "Check if VTK should be generated." << std::endl;
			std::stringstream ss;
			ss << cfg.output_folder << cfg.identifier << "_" << std::setw(4) << std::setfill('0') << std::to_string(vtkcount + 1) << '_' << std::to_string((size_t)timestep) << ".vtk";
			std::string vtk_path = ss.str();
			vtk::to_vtk(ss.str().c_str(), cube);
			std::string json_path; // json path

			if (cfg.log_transitions) cube->flush_log_file();

			// write to analysis file
			if (cfg.generate_analysis_files)
			{
				std::cout << "Beginning analysis..." << std::endl;
				analyze.load_lattice(cube);
				ss.str(std::string());
				//ss << cfg.output_folder << cfg.identifier << "_" << std::setw(4) << std::setfill('0') << std::to_string(vtkcount + 1) << '_' << std::to_string((size_t)timestep) << "_analysis.txt";
				ss << cfg.output_folder << cfg.identifier << "_" << std::setw(4) << std::setfill('0') << std::to_string(vtkcount + 1) << '_' << std::to_string((size_t)timestep) << "_analysis.json";
				
				// josn file
				json_path = ss.str();
				analyze.save_analysis_to_json(ss.str().c_str(), timestep);
			}
			bundle_and_cleanup(vtk_path, json_path);
			small_zips.push_back(vtk_path.substr(0, vtk_path.size() - 4) + ".zip");

			++vtkcount;
			++curr_checkpoint;

			if (cfg.max_timestep <= 0 && curr_checkpoint >= checkpoints.size()) break;
		}
		else if (cfg.checkpoint_interval > 0 && timestep >= next_checkpoint)
		{
			std::stringstream ss;
			ss << cfg.output_folder << cfg.identifier << "_" << std::setw(4) << std::setfill('0') << std::to_string(vtkcount + 1) << '_' << std::to_string((size_t)timestep) << ".vtk";
			std::string vtk_path = ss.str();
			vtk::to_vtk(vtk_path.c_str(), cube);

			std::string json_path;
			// // Add VTK to zip
			// zip_source_t* src = zip_source_file(archive, vtk_path.c_str(), 0, -1);
			// if (src) {
			// 	auto idx = zip_file_add(archive, std::filesystem::path(vtk_path).filename().c_str(), src, ZIP_FL_OVERWRITE);
			// 	zip_set_file_compression(archive, idx, ZIP_CM_DEFLATE, 9);
			// }

			if (cfg.log_transitions) cube->flush_log_file();

			if (cfg.generate_analysis_files)
			{
				std::cout << "Grain count: " << analyze.get_grain_count() << std::endl;
				std::cout << "Beginning analysis..." << std::endl;
				analyze.load_lattice(cube);

				ss.str(std::string());
				ss << cfg.output_folder << cfg.identifier << "_" << std::setw(4) << std::setfill('0') << std::to_string(vtkcount + 1) << '_' << std::to_string((size_t)timestep) << "_analysis.json";
				
				json_path = ss.str();
				analyze.save_analysis_to_json(json_path.c_str(), timestep);

				// // Add JSON to zip
				// zip_source_t* jsrc = zip_source_file(archive, json_path.c_str(), 0, -1);
				// if (jsrc) {
				// 	auto idx = zip_file_add(archive, std::filesystem::path(json_path).filename().c_str(), jsrc, ZIP_FL_OVERWRITE);
				// 	zip_set_file_compression(archive, idx, ZIP_CM_DEFLATE, 9);
				// }
			}
			bundle_and_cleanup(vtk_path, json_path);
			small_zips.push_back(vtk_path.substr(0, vtk_path.size() - 4) + ".zip");

			++vtkcount;
			next_checkpoint += cfg.checkpoint_interval;
		}

		// Break if the max timestep is reached.
		if (cfg.max_timestep > 0 && timestep >= cfg.max_timestep) break;

		// break if grain size is less than 5
		if (cfg.generate_analysis_files && analyze.get_grain_count() > 0 && analyze.get_grain_count() <= 3905) {
			std::cout << "Grain count reached " << analyze.get_grain_count() << ", stopping." << std::endl;
			break;
		}

		//break; // testing 
		// break if there are 5 grains left
		// count how many grains are rleft

	}

	if (cfg.log_transitions) cube->stop_logging_transitions();

	// Close zip after the loop — this is when everything gets written
	// if (zip_close(archive) < 0){
	// 	std::cerr << "zip_close failed: " << zip_strerror(archive) << std::endl;
	// }
	// merge into large zip

	std::cout << "Small zips to merge: " << small_zips.size() << std::endl;
	for (const auto& z : small_zips){
		std::cout << "  " << z << " exists=" << std::filesystem::exists(z) << std::endl;
	}
	std::string final_zip = cfg.output_folder + cfg.identifier + "_all.zip";
	merge_zips(small_zips, final_zip);



	// Commented out code for octree testing

	/*lattice_t cube(1);
	size_t size = 100;
	unsigned char height = 6;
	voxel_t *vlist = new voxel_t[size * size * size];
	octree3_t tree(128, log2(128));
	for (size_t z = 0; z < size; ++z)
	{
		for (size_t y = 0; y < size; ++y)
			for (size_t x = 0; x < size; ++x) {
				tree.delta(x, y, z, 1);
				vlist[x + y * size + z * size * size].activity = 1;
			}
	}

	for (size_t i = 0; i < 10000; ++i)
	{
		coord_t x, y, z;
		x = cube.rng(0, size);
		y = cube.rng(0, size);
		z = cube.rng(0, size);

		std::cout << "=================" << std::endl;
		std::cout << "real: " << x << ", " << y << ", " << z << std::endl;

		tree.delta(x, y, z, 1);
		vlist[x + y * size + z * size * size].activity = 1;
		std::cout << "setting vindex " << x + y * size + z * size * size << " to 1 " << std::endl;

		coord_t vx, vy, vz;
		tree.get_voxel_from_sum_activity(&vx, &vy, &vz, 1, vlist, 100);

		std::cout << "pred: " << vx << ", " << vy << ", " << vz << std::endl;

		if (x != vx || y != vy || z != vz)
		{
			std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!! MISMATCH!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		}

		tree.delta(x, y, z, -1);
		vlist[x + y * size + z * size * size].activity = 0;
	}*/

	//tree.dump();
	//tree.dump_level(height - 1);

	/*coord_t vx, vy, vz;
	tree.get_voxel_from_sum_activity(&vx, &vy, &vz, size * size * size / 2, vlist);
	std::cout << vx << ", " << vy << ", " << vz << std::endl;*/
}