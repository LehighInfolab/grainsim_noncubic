#pragma once

// Convert the 1d array of spin into 3d structure for mat format export

#include <string>
#include <vector>
#include <fstream>
#include <cstdint>
#include <iostream>
#include <algorithm>

#include <matio.h>

#include "lattice.h"
#include "types.h"

// for actual grain writing and the mobility cube
template <typename T> struct mat_type; 
template <> struct mat_type<uint8_t>
{
    static constexpr enum matio_classes cls = MAT_C_UINT8;
    static constexpr enum matio_types   typ = MAT_T_UINT8;
};
template <> struct mat_type<uint32_t>
{
    static constexpr enum matio_classes cls = MAT_C_UINT32;
    static constexpr enum matio_types   typ = MAT_T_UINT32;
};

// template for 3d cubes 
template <typename T>
// 3d representation of array 
// index = x + y*nx + z*nx*ny
struct cube3d_t
{
    coord_t nx = 0, ny =0, nz = 0;
    std::vector<T> data;

    void alloc(const lattice_t *lat)
    {
        nx = lat->side_length_x;
        ny = lat->side_length_y;
        nz = lat->side_length_z;
        data.resize((size_t)nx * ny * nz);
    }

    size_t index(coord_t x, coord_t y, coord_t z) const
    {
        return x + (size_t)y * nx + (size_t)z * nx * ny;
    }

    size_t nbytes() const { return data.size() * sizeof(T); }
};

class array3d
{
    static_assert(sizeof(spin_t) == 4, "array3d assumes 4-byte spin_t");

    // read only boundary transition check
    static bool boundary_transformed(const lattice_t *lattice, spin_t a, spin_t b)
    {
        spin_t lo = std::min(a, b), hi = std::max(a, b);
        auto &bmap = lattice->boundary_tracker.boundary_map;
        auto sm = bmap.find(lo);
        if (sm == bmap.end()) {return false;}
        auto lg = sm->second.find(hi);
        if (lg == sm->second.end()) {return false;}
		return lg->second->transformed;
    }

public:
    

    // convert 1d array to 3d shape lattice
    static cube3d_t<spin_t> to_cube(lattice_t *lattice)
    {
        // cube3d_t cube;
        // cube.nx = lattice->side_length_x;
        // cube.ny = lattice->side_length_y;
        // cube.nz = lattice->side_length_z;
        // cube.data.resize((size_t)cube.nx * cube.ny * cube.nz);
        cube3d_t<spin_t> cube;
        cube.alloc(lattice);

        // for (coord_t z = 0; z < cube.nz; ++z) 
        //     for (coord_t y = 0; y < cube.ny; ++y)
        //         for (coord_t x = 0; x < cube.nx; ++x)
        //         {
        //             size_t idx = x + (size_t)y * cube.nx + (size_t)z * cube.nx * cube.ny;
        //             cube.data[idx] = lattice->voxels[idx].spin;
        //         }

        const size_t n = cube.data.size();
        for (size_t i = 0; i < n; ++i)
        {
            cube.data[i] = lattice->voxels[i].spin;
        }
        return cube;
    }

    // to mobility cube 
    // 0 -> interior voxels, all 26 neighbors share this 
    // 1 -> untransformed boundary, touch another grain but not transformed 
    // 2 -> transformed bounary
    static cube3d_t<uint8_t> to_cube_mobility(lattice_t *lattice)
    {
        cube3d_t<uint8_t> cube;
        cube.alloc(lattice);

        for (coord_t z = 0; z < cube.nz; ++z)
			for (coord_t y = 0; y < cube.ny; ++y)
				for (coord_t x = 0; x < cube.nx; ++x)
				{
					spin_t s = lattice->voxel_at(x, y, z)->spin;
					uint8_t state = 0; // interior as default
					for (int n = 0; n < NEIGH_COUNT; ++n)
					{
						spin_t ns = lattice->neighbor_at(x, y, z, n)->spin;
						if (ns == s) continue;          
						state = 1;        // on untransformed oundary
						if (boundary_transformed(lattice, s, ns))
						{
							state = 2;                   // transformed boundary
							break;
						}
					}
					// size_t idx = x + (size_t)y * cube.nx + (size_t)z * cube.nx * cube.ny;
					// cube.data[idx] = state;
                    cube.data[cube.index(x, y, z)] = state;
				}
        
        return cube;
    }

    // ===========================================================
    // convert a cube to .mat format
    // ===========================================================

    template <typename T>
    static void to_mat(const char *fname, const cube3d_t<T> &cube,  const std::string &var_name = "grains")
    {
        std::cout << "Writing MAT array to " << fname << std::endl;

        const size_t MAT5_LIMIT = (size_t)2 * 1024 * 1024 * 1024 - 16 * 1024;
        enum mat_ft version = (cube.nbytes() >= MAT5_LIMIT) ? MAT_FT_MAT73: MAT_FT_MAT5;

        mat_t *matfp = Mat_CreateVer(fname, nullptr, version);
        if (matfp == nullptr)
        {
            std::cerr << "array3d: Mat_CreateVer failed for \"" << fname
                      << "\" (is matio built with MAT73 support?)." << std::endl;
            return;
        }

        size_t dims[3] = {(size_t)cube.nx, (size_t)cube.ny, (size_t)cube.nz};

        matvar_t *matvar = Mat_VarCreate(
            var_name.c_str(), mat_type<T>::cls, mat_type<T>::typ, 3, dims,
            const_cast<T *>(cube.data.data()), 0);
            
        if (matvar == nullptr)
        {
            std::cerr << "array3d: Mat_VarCreate failed." << std::endl;
            Mat_Close(matfp);
            return;
        }

        if (Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB) != 0)
            std::cerr << "array3d: Mat_VarWrite failed." << std::endl;

        Mat_VarFree(matvar);
        Mat_Close(matfp);

        std::cout << "Done writing .mat!" << std::endl;
    }

    // static void to_mat(const char *fname, const cube3d_t &cube, const std::string &var_name = "grains")
    // {
    //     std::cout << "Writing MAT array to " << fname << std::endl;
    //     size_t nbytes = cube.data.size() * sizeof(spin_t);

    //     // if smaller than 2gb 
    //     const size_t MAT5_LIMIT = (size_t)2 * 1024 * 1024 * 1024 - (16 * 1024);
    //     enum mat_ft version = (nbytes >= MAT5_LIMIT) ? MAT_FT_MAT73 : MAT_FT_MAT5;

    //     mat_t *matfp = Mat_CreateVer(fname, nullptr, version);
	// 	if (matfp == nullptr)
	// 	{
	// 		std::cerr << "array3d: Mat_CreateVer failed for \"" << fname
	// 				  << "\" (is matio built with MAT73 support?)." << std::endl;
	// 		return;
	// 	}


    //     // column major for mat 
    //     size_t dims[3] = {(size_t)cube.nx, (size_t)cube.ny, (size_t)cube.nz};


    //     matvar_t *matvar = Mat_VarCreate(
	// 		var_name.c_str(), MAT_C_UINT32, MAT_T_UINT32, 3, dims,
	// 		const_cast<spin_t *>(cube.data.data()), 0);
            
	// 	if (matvar == nullptr)
	// 	{
	// 		std::cerr << "array3d: Mat_VarCreate failed." << std::endl;
	// 		Mat_Close(matfp);
	// 		return;
	// 	}

	// 	if (Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB) != 0)
	// 	{
	// 		std::cerr << "array3d: Mat_VarWrite failed." << std::endl;
	// 	}

	// 	Mat_VarFree(matvar);
	// 	Mat_Close(matfp);

	// 	std::cout << "Done writing .mat!" << std::endl;
    // }

};