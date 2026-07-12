#pragma once

// Convert the 1d array of spin into 3d structure for mat format export

#include <string>
#include <vector>
#include <fstream>
#include <cstdint>
#include <iostream>

#include <matio.h>

#include "lattice.h"
#include "types.h"

class array3d
{
    static_assert(sizeof(spin_t) == 4, "array3d assumes 4-byte spin_t");
public:
    // 3d representation of array 
    // index = x + y*nx + z*nz*nz
    struct cube3d_t
    {
        coord_t nx = 0, ny =0, nz = 0;
        std::vector<spin_t> data;
    };

    static cube3d_t to_cube(lattice_t *lattice)
    {
        cube3d_t cube;
        cube.nx = lattice->side_length_x;
        cube.ny = lattice->side_length_y;
        cube.nz = lattice->side_length_z;
        cube.data.resize((size_t)cube.nx * cube.ny * cube.nz);

        for (coord_t z = 0; z < cube.nz; ++z) 
            for (coord_t y = 0; y < cube.ny; ++y)
                for (coord_t x = 0; x < cube.nx; ++x)
                {
                    size_t idx = x + (size_t)y * cube.nx + (size_t)z * cube.nx * cube.ny;
                    cube.data[idx] = lattice->voxels[idx].spin;
                }

        return cube;
    }

    // convert a cube to .mat format
    static void to_mat(const char *fname, const cube3d_t &cube, const std::string &var_name = "grains")
    {
        std::cout << "Writing MAT array to " << fname << std::endl;
        size_t nbytes = cube.data.size() * sizeof(spin_t);

        // if smaller than 2gb 
        const size_t MAT5_LIMIT = (size_t)2 * 1024 * 1024 * 1024 - (16 * 1024);
        enum mat_ft version = (nbytes >= MAT5_LIMIT) ? MAT_FT_MAT73 : MAT_FT_MAT5;

        mat_t *matfp = Mat_CreateVer(fname, nullptr, version);
		if (matfp == nullptr)
		{
			std::cerr << "array3d: Mat_CreateVer failed for \"" << fname
					  << "\" (is matio built with MAT73 support?)." << std::endl;
			return;
		}


        // column major for mat 
        size_t dims[3] = {(size_t)cube.nx, (size_t)cube.ny, (size_t)cube.nz};


        matvar_t *matvar = Mat_VarCreate(
			var_name.c_str(), MAT_C_UINT32, MAT_T_UINT32, 3, dims,
			const_cast<spin_t *>(cube.data.data()), 0);
            
		if (matvar == nullptr)
		{
			std::cerr << "array3d: Mat_VarCreate failed." << std::endl;
			Mat_Close(matfp);
			return;
		}

		if (Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB) != 0)
		{
			std::cerr << "array3d: Mat_VarWrite failed." << std::endl;
		}

		Mat_VarFree(matvar);
		Mat_Close(matfp);

		std::cout << "Done writing .mat!" << std::endl;
    }

};