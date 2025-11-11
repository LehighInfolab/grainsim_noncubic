
#include <iostream>
#include <vector>
#include <sstream>
#include <iomanip>

#include "vtk.h"
#include "lattice.h"
#include "octree3.h"
#include "debug_timer.h"
#include "config.h"
#include "analysis.h"

void print3DArray(const std::vector<int>& arr, int dim_x, int dim_y, int dim_z) {
    if (arr.size() != static_cast<size_t>(dim_x * dim_y * dim_z)) {
        std::cerr << "Error: array size does not match dimensions!\n";
        return;
    }

    for (int z = 0; z < dim_z; ++z) {
        std::cout << "z = " << z << ":\n";
        for (int y = 0; y < dim_y; ++y) {
            for (int x = 0; x < dim_x; ++x) {
                int index = z * (dim_x * dim_y) + y * dim_x + x;
                std::cout << arr[index] << " ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }
}

// g++ -O3 grainsim/tests.cpp -o tests.out -static

int main(int argc, char *argv[]) {
    // Load the config file.
	config_t cfg;
	cfg.load_config();

    // Create a lattice for testing
    

    coord_t dim_x = 36;
    coord_t dim_y= 36; 
    coord_t dim_z = 128;

    lattice_t noncube(dim_x, dim_y, dim_z);

    voxel_t *vlist = new voxel_t[dim_x* dim_y*dim_z];
    std::vector<int> arr(dim_x* dim_y*dim_z);
    octree3_t tree(dim_z, log2(dim_z));
    
    for (size_t z = 0; z < dim_z; ++z) {
        for (size_t y = 0; y < dim_y; ++y) {
            for (size_t x = 0; x < dim_x; ++x) {
                tree.delta(x, y, z, 1);
				vlist[x + y * dim_y + z * dim_z * dim_z].activity = 1;
                arr[x + y * dim_y + z * dim_z * dim_z] = 1;
            }
        }
    }

    print3DArray(arr, dim_x, dim_y, dim_z);
    return 0;

}


