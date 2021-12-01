//
// Created by johannes on 12/1/21.
//

#ifndef POWER_SPECTRUM_TRANSFORMATIONS_H
#define POWER_SPECTRUM_TRANSFORMATIONS_H

#include "aweights.h"
#include "blitz/array.h"

/**
 * Wrap the point i on a grid with size n_grid using if else statements.
 *
 * @param i Coordinate to wrap.
 * @param n_grid Wrap modulo n_grid.
 */
int wrap_if_else(int i, int n_grid) {
    if (i >= n_grid)
        return i - n_grid;
    if (i < 0)
        return i + n_grid;
    return i;
}


/**
 * Wrap the point i on a grid with size n_grid using the % operator.
 *
 * @param i Coordinate to wrap.
 * @param n_grid Wrap modulo n_grid.
 */
int wrap_modulo(int i, int n_grid) { return (i + n_grid) % n_grid; }
#endif // POWER_SPECTRUM_TRANSFORMATIONS_H


/**
 * Return the grid coordinate for a euclidean coordinate in the interval
 * [-0.5, 0.5].
 *
 * @param coord: Coordinate.
 * @param n_grid: Size of the grid.
 */
template <typename real_t> real_t grid_coordinate(real_t coord, int n_grid) {
    return n_grid * (coord + 0.5);
}


/**
 * Return the rank that a particle belongs to in the mass grid.
 *
 * @param grid_coordinate: X-coordinate of the particle in the grid.
 * @param starting_indices: Pointer to array that contains the starting indices
 * of all MPI ranks.
 * @param size_starting_indices: Size of starting_indices array.
 *
 */
template <int Order, typename real_t>
int particle_rank(real_t grid_coordinate, std::ptrdiff_t *starting_indices,
                  int size_starting_indices) {
    auto w = AssignmentWeights<Order>(grid_coordinate);
    int rank_index = 0;
    for (auto i = 0; i < size_starting_indices; i++) {
        if (w.i < starting_indices[i]) {
            rank_index = i - 1;
            break;
        }
    }
    return wrap_if_else(rank_index, size_starting_indices);
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
