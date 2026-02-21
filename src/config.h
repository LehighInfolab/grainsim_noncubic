#pragma once

#include <vector>
#include <unistd.h>
#include <sstream>
#include <fstream>

#include "types.h"

struct config_t
{
	std::string initial_state_path;
	std::string identifier;
	std::string output_folder;
	std::string checkpoints;
	double checkpoint_interval = -1;
	double max_timestep = -1;
	activ_t default_mobility = 0.002, transitioned_mobility = 0.04;
	double transition_interval = 0;
	size_t transition_count = 0;
	double scale_multiplier = 1;
	double propagation_chance = 0.95;
	bool use_potential_energy = false;
	int const_grain_count = 0;
	bool log_transitions = true;
	double propagation_ratio = 0;
	bool generate_analysis_files = true;
	// Z axis for the transition, by default it's 0
	int z_propagation_plane = 0; 

	// increase in transition time count in interval
	size_t transition_time_count_increase = 0;
	// switch to turn rampup on or off
	bool ramp_up_enabled = false;
	// high nucleation rate limit
	double high_nucleation_rate_limit = 0;

	void checkpoints_to_vector(std::vector<double> *checkpoint_vector)
	{
		std::istringstream ss(checkpoints);
		std::string word;
		while (ss >> word)
		{
			checkpoint_vector->push_back(std::stoul(word));
		}
	}
	// Modify for slurm batch process with input and output override
	void load_config(const std::string &config_path, 
		const std::string &initial_state_override, 
		const std::string &output_override,
		int transition_count_override)
	{
		//std::ifstream cfgfile("grainsim_config.txt");
		std::ifstream cfgfile(config_path);
		if (!cfgfile) {
			std::cerr << "Error: could not open config file \"" << config_path << "\"\n";
			std::exit(1);
		}
		
		std::string line;
		while (std::getline(cfgfile, line))
		{
			std::istringstream ss(line);
			std::string key, equals, value, tempv;

			if(!(ss >> key)) continue;
			if (key[0] == '#') continue;

			if (!(ss >> equals)) continue;
			if (!(ss >> value)) continue;
			while (ss >> tempv)
			{
				value.append(" ");
				value.append(tempv);
			}

			if (key == "INITIAL_STATE_FILE")
			{
				initial_state_path = value;
			}
			else if (key == "OUTPUT_FOLDER")
			{
				output_folder = value;
			}
			else if (key == "IDENTIFIER")
			{
				identifier = value;
			}
			else if (key == "CHECKPOINTS")
			{
				checkpoints = value;
			}
			else if (key == "PERIODIC_CHECKPOINT_INTERVAL")
			{
				checkpoint_interval = std::stod(value);
			}
			else if (key == "MAX_TIMESTEP")
			{
				max_timestep = std::stod(value);
			}
			else if (key == "DEFAULT_MOBILITY")
			{
				default_mobility = std::stod(value);
			}
			else if (key == "TRANSITIONED_MOBILITY")
			{
				transitioned_mobility = std::stod(value);
			}
			else if (key == "TRANSITION_INTERVAL")
			{
				transition_interval = std::stod(value);
			}
			else if (key == "TRANSITION_COUNT")
			{
				transition_count = std::stoi(value);
			}
			else if (key == "PROPAGATION_CHANCE")
			{
				propagation_chance = std::stod(value);
			}
			else if (key == "USE_POTENTIAL_ENERGY")
			{
				use_potential_energy = value == "true";
			}
			else if (key == "SCALE_MULTIPLIER")
			{
				scale_multiplier = std::stod(value);
			}
			else if (key == "LOG_BOUNDARY_TRANSITIONS")
			{
				log_transitions = value == "true";
			}
			else if (key == "CONST_GRAIN_COUNT")
			{
				const_grain_count = std::stoi(value);
			}
			else if (key == "PROPAGATION_RATIO")
			{
				propagation_ratio = std::stod(value);
			}
			else if (key == "GENERATE_ANALYSIS_FILES")
			{
				generate_analysis_files = value == "true";
			}
			// add config: z axis
			else if (key == "Z_PROP_PLANE"){
				z_propagation_plane = std::stoi(value);
				std::cout << "[INFO] Z_PROP_PLANE: " << z_propagation_plane << std::endl;
			}
			else if (key == "TRANSITION_TIME_COUNT_INCREASE")
			{
				transition_time_count_increase = std::stoul(value);
			}
			else if (key == "RAMP_UP_ENABLED")
			{
				ramp_up_enabled = value == "true";
			}
			else if (key == "HIGH_NUCLEATION_RATE_LIMIT")
			{
				high_nucleation_rate_limit = std::stod(value);
			}
			else
			{
				std::cout << "Warning: Unknown config key \"" << key << "\"." << std::endl;
			}
		}

		// Override with function config
		if (!initial_state_override.empty()) {
			initial_state_path = initial_state_override;
		}
		if (!output_override.empty()) {
			output_folder = output_override;
		}
		// overrite transition count
		if (transition_count_override>= 0) {
			transition_count = transition_count_override;
		}
		// make sure ends with slash
		if (!output_folder.empty()) {
			char last = output_folder.back();
			if (last != '/' && last != '\\'){
				output_folder.push_back('/');
			}
		}
	}
};