//
// Created by johannes on 12/1/21.
//

#ifndef POWER_SPECTRUM_TRANSFORMATIONS_H
#define POWER_SPECTRUM_TRANSFORMATIONS_H

#include "aweights.h"
#include "blitz/array.h"

/*
 * Return the grid coordinate from a euclidean coordinate.
 *
 * @param coord coordinate.
 * @param n_grid size of the grid.
 */
template <typename real_t> real_t grid_coordinate(real_t coord, int n_grid) {
    return n_grid * (coord + 0.5);
}

/**
 * Return the rank that a particle belongs to in the mass grid.
 *
 * @param x_coordinate: X-coordinate of the particle.
 * @param n_grid: Number of grid cells in the x-direction.
 * @param starting_indices: Pointer to array that contains the starting indices
 * of all MPI ranks.
 * @param size_starting_indices: Size of starting_indices array.
 *
 */
template <int Order, typename real_t>
int particle_rank(real_t x_coordinate, int n_grid,
                  std::ptrdiff_t *starting_indices, int size_starting_indices) {
    auto w = AssignmentWeights<Order>(grid_coordinate(x_coordinate, n_grid));
    int rank_index = 0;
    for (auto i = 0; i < size_starting_indices; i++) {
        if (w.i < starting_indices[i]) {
            rank_index = i - 1;
            break;
        }
    }
    return rank_index;
}

template <typename real_t>
blitz::Array<real_t, 2>
project_3d_to_2d(const blitz::Array<real_t, 3> &grid_3d) {
    int n_grid = grid_3d.extent(blitz::firstDim);
    blitz::thirdIndex k;
    blitz::Array<real_t, 2> grid_2d(n_grid, n_grid);

    grid_2d = max(grid_3d, k);
    return grid_2d;
}
#endif // POWER_SPECTRUM_TRANSFORMATIONS_H
