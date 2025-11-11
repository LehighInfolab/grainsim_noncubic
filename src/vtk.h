#pragma once

#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>

#include "lattice.h"
#include "types.h"

// A simple class that allows for reading from and writing to .vtk files.
class vtk
{
private:
	// Determine if a string starts with another string.
	static bool str_starts_with(std::string str, const char *exp)
	{
		size_t index = 0;
		while (exp[index] != 0)
		{
			if (str[index] != exp[index]) return false;
			++index;
		}
		return true;
	}

	// Determine if a string ends with another string.
	static bool str_ends_with(std::string str, const char *exp)
	{
		size_t explen = strlen(exp);

		for (size_t i = 0; i < explen; ++i)
		{
			if (str[str.length() - explen + i] != exp[i]) return false;
		}
		return true;
	}

public:
	// Create a lattice object from a .vtk file.
	static lattice_t *from_vtk(const char *fname, bool init=true)
	{
		std::cout << "Loading VTK file " << fname << std::endl;

		lattice_t *new_cube;

		std::ifstream vtkfile(fname);
		std::string line;
		spin_t index = 0;
		char load_state = 0;
		bool end_loop = false;
		while (std::getline(vtkfile, line))
		{
			// Cases govern status of parser.
			// Read cube dimension sizes, then skip to CELL_DATA section, then read voxel spins.
			switch (load_state)
			{
			case 0:
				if (str_starts_with(line, "DIMENSIONS"))
				{
					std::istringstream ss(line);
					std::string word;
					ss >> word; // skip "DIMENSIONS"
					ss >> word;
					//vtk header format: DIMENSIONS 151 151 151 
					coord_t dim_x = std::stoi(word) -1;
					ss >> word;
					coord_t dim_y = std::stoi(word) -1;
					ss >> word;
					coord_t dim_z = std::stoi(word) -1;
					std::cout << "vtk x, y, z: " << dim_x << ", " << dim_y << ", " << dim_z << std::endl;
					//new_cube = new lattice_t(dim);
					// updated constructor
					new_cube = new lattice_t(dim_x, dim_y, dim_z);
					++load_state;
				}
				break;
			case 1:
				if (str_starts_with(line, "CELL_DATA")) ++load_state;
				break;
			case 2:
				if (line[0] >= '0' && line[0] <= '9')
				{
					new_cube->voxels[index++].spin = std::stoul(line);
					++load_state;
				}
				break;
			case 3:
				if (line[0] >= '0' && line[0] <= '9')
				{
					new_cube->voxels[index++].spin = std::stoul(line);
				}
				else end_loop = true;
				break;
			}

			if (end_loop) break;
		}

		vtkfile.close();

		std::cout << "Done loading!" << std::endl;
		if (init)
		{
			new_cube->init();
		}

		return new_cube;
	}

	// Save a lattice object to a .vtk file.
	static void to_vtk(const char *fname, lattice_t *lattice)
	{
		std::ofstream vtkfile(fname);

		std::cout << "Writing to " << fname << std::endl;

		vtkfile << "# vtk DataFile Version 2.0\n data set from May6 1\nASCII\nDATASET RECTILINEAR_GRID\n";
		vtkfile << "DIMENSIONS " << (lattice->side_length_x + 1) << " " << (lattice->side_length_y + 1) << " " << (lattice->side_length_z + 1) << " \n";

		vtkfile << "X_COORDINATES " << (lattice->side_length_x + 1) << " Float \n";
		for (size_t i = 0; i < lattice->side_length_x + 1; ++i)
		{
			vtkfile << i << '\n';
		}
		vtkfile << "Y_COORDINATES " << (lattice->side_length_y + 1) << " Float \n";
		for (size_t i = 0; i < lattice->side_length_y + 1; ++i)
		{
			vtkfile << i << '\n';
		}
		vtkfile << "Z_COORDINATES " << (lattice->side_length_z + 1) << " Float \n";
		for (size_t i = 0; i < lattice->side_length_z + 1; ++i)
		{
			vtkfile << i << '\n';
		}
		vtkfile << "CELL_DATA " << (lattice->side_length_x * lattice->side_length_y * lattice->side_length_z) << " \n";
		vtkfile << "SCALARS GrainIDs int  1\nLOOKUP_TABLE default\n";
		for (size_t i = 0; i < (lattice->side_length_x * lattice->side_length_y * lattice->side_length_z); ++i)
		{
			vtkfile << lattice->voxels[i].spin << '\n';
		}
		
		vtkfile.close();
	}

	// Create a lattice object from a .ph file.
	static lattice_t *from_ph(const char *fname, bool init=true)
	{
		std::cout << "Loading PH file " << fname << std::endl;

		lattice_t *new_cube;

		std::ifstream phfile(fname);
		std::string line;
		spin_t index = 0;
		char load_state = 0;
		bool end_loop = false;
		while (std::getline(phfile, line))
		{
			// Cases govern status of parser.
			switch (load_state)
			{
			case 0:
			{
				std::istringstream ss(line);
				std::string word;
				ss >> word;
				//ph header format:      256     256     1024
				coord_t dim_x = std::stoi(word);
				ss >> word;
				coord_t dim_y = std::stoi(word);
				ss >> word;
				coord_t dim_z = std::stoi(word);
				std::cout << "x, y, z: " << dim_x << ", " << dim_y << ", " << dim_z << std::endl;
				//new_cube = new lattice_t(dim);
				// updated constructor
				new_cube = new lattice_t(dim_x, dim_y, dim_z);
				++load_state;

				break;
			}
			case 1:
			case 2:
				++load_state;
				break;
			case 3:
				if (line[0] >= '0' && line[0] <= '9')
				{
					new_cube->voxels[index++].spin = std::stoul(line);
				}
				else end_loop = true;
				break;
			}

			if (end_loop) break;
		}

		phfile.close();

		std::cout << "Done loading!" << std::endl;
		if (init)
		{
			new_cube->init();
		}

		return new_cube;
	}

	// Load a file and autodetect the correct load function to use.
	static lattice_t *from_file(const char *fname, bool init=true)
	{
		std::string str = std::string(fname);

		if (str_ends_with(str, ".vtk"))
		{
			return from_vtk(fname, init);
		}
		else if (str_ends_with(str, ".ph"))
		{
			return from_ph(fname, init);
		}
		else
		{
			std::cout << "Error: Unrecognized file format." << std::endl;
			exit(0);
		}
	}

	static lattice_t *scale_lattice(lattice_t *lat, double multiplier, bool init=true)
	{
		std::cout << "Scaling lattice..." << std::endl;

		lattice_t *new_cube = new lattice_t(lat->side_length_x * multiplier, lat->side_length_y * multiplier, lat->side_length_z * multiplier);

		for (int z = 0; z < new_cube->side_length_z; ++z)
			for (int y = 0; y < new_cube->side_length_y; ++y)
				for (int x = 0; x < new_cube->side_length_x; ++x)
				{
					new_cube->voxels[
							(x) +
							(y) * new_cube->side_length_y +
							(z) * new_cube->side_length_z * new_cube->side_length_z].spin
						= lat->voxels[
							(int)(x / multiplier) +
							(int)(y / multiplier) * lat->side_length_y +
							(int)(z / multiplier) * lat->side_length_z * lat->side_length_z].spin;
				}

		if (init)
		{
			new_cube->init();
		}

		return new_cube;
	}
};