#pragma once

// The number of neighbors a voxel can have.
const char NEIGH_COUNT = 26;

// The data type that holds a grain's spin (ID).
typedef unsigned int spin_t;
// The data type that holds an activity (probability).
typedef double activ_t;
// The data type that holds a single dimension within the lattice (must be signed to allow for wrapping).
typedef int coord_t;

// selection type for sweeping mechanism
enum sweep_selection_t
{
    SWEEP_SMALLEST = 0,
    SWEEP_LARGEST = 1,
    SWEEP_RANDOM = 2
};