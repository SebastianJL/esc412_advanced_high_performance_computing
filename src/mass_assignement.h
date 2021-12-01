//
// Created by johannes on 11/12/21.
//

#ifndef PROJECT_MASS_ASSIGNEMENT_H
#define PROJECT_MASS_ASSIGNEMENT_H

#include "aweights.h"
#include "blitz/array.h"
#include "io_utils.h"
#include "transformations.h"

/*
 * Return the nearest grid point for a single coordinate.
 *
 * @param coord coordinate.
 * @param n_grid size of the grid.
 */
template<typename real_t>
real_t ngp(real_t coord, int n_grid) {
    return std::floor(n_grid * (coord + 0.5));
}
/*
 * Return a grid of mass densities given a particle distribution.
 *
 * The algorithm used is nearest grid point.
 *
 * @param particles Particles to distribute
 * @param n_grid Grid size in 3d.
 */
template<typename real_t>
blitz::Array<real_t, 3> assign_mass_ngp(blitz::Array<real_t, 2> particles,
                                        int n_grid) {
    blitz::Array<real_t, 3> res(n_grid, n_grid, n_grid);
    for (auto i = 0; i < particles.extent(blitz::firstDim); ++i) {
        blitz::TinyVector<real_t, 3> grid_point;
        for (auto j = 0; j < particles.extent(blitz::secondDim); ++j) {
            real_t p = particles(i, j);
            grid_point(j) = ngp(p, n_grid);
        }
        real_t mass = 1.;
        res(grid_point) += mass;
    }
    return res;
}
/*
 * Return a grid of mass densities given a particle distribution.
 *
 * @param particles Particles to distribute.
 * @param n_grid Grid size in 3d.
 * @param wrap Wrapping function to map values from [-n_grid, 2*n_grid) to the
 * range [0, n_grid). It should be equivalent to (i + n_grid) % n_grid.
 * @param out Grid for assigning the masses.
 * @param order Order of the mass assignment weights.
 */
template <typename real_t, int Order>
void assign_mass(blitz::Array<real_t, 2> particles, int n_grid,
                 int (*wrap)(int, int), blitz::Array<real_t, 3> out) {
    out = 0;
    for (auto row = 0; row < particles.extent(blitz::firstDim); ++row) {
        real_t p = particles(row, 0);
        real_t p_grid = grid_coordinate(p, n_grid);
        AssignmentWeights<Order> wx(p_grid);

        p = particles(row, 1);
        p_grid = grid_coordinate(p, n_grid);
        AssignmentWeights<Order> wy(p_grid);

        p = particles(row, 2);
        p_grid = grid_coordinate(p, n_grid);
        AssignmentWeights<Order> wz(p_grid);

        for (auto i = 0; i < Order; ++i) {
            for (auto j = 0; j < Order; ++j) {
                for (auto k = 0; k < Order; ++k) {
                    int x = wrap(wx.i + i, n_grid);
                    int y = wrap(wy.i + j, n_grid);
                    int z = wrap(wz.i + k, n_grid);
                    out(x, y, z) += wx.H[i] * wy.H[j] * wz.H[k];
                }
            }
        }
    }
}
/*
 * Return a grid of mass densities given a particle distribution.
 *
 * Use a grid bigger with margins instead of using a modulo operation.
 *
 * @param particles Particles to distribute
 * @param n_grid Grid size in 3d.
 * @param order Order of the mass assignment weights.
 */
template <typename real_t, int Order>
blitz::Array<real_t, 3>
assign_mass_with_margins(blitz::Array<real_t, 2> particles, int n_grid) {
    using blitz::Array;
    using blitz::firstDim;
    using blitz::Range;

    int margin = (Order - 1);
    int n_grid_m = n_grid + 2 * margin;
    blitz::Array<real_t, 3> res_m(n_grid_m, n_grid_m, n_grid_m);
    res_m = 0;

    blitz::Range inside(blitz::Range(margin, n_grid + margin - 1));
    blitz::Array<real_t, 3> res = res_m(inside, inside, inside);

    for (auto row = 0; row < particles.extent(blitz::firstDim); ++row) {
        real_t p = particles(row, 0);
        real_t p_grid = grid_coordinate(p, n_grid);
        AssignmentWeights<Order> wx(p_grid);

        p = particles(row, 1);
        p_grid = grid_coordinate(p, n_grid);
        AssignmentWeights<Order> wy(p_grid);

        p = particles(row, 2);
        p_grid = grid_coordinate(p, n_grid);
        AssignmentWeights<Order> wz(p_grid);

        for (auto i = 0; i < Order; ++i) {
            for (auto j = 0; j < Order; ++j) {
                for (auto k = 0; k < Order; ++k) {
                    int x = wx.i + i + margin;
                    int y = wy.i + j + margin;
                    int z = wz.i + k + margin;
                    res_m(x, y, z) += wx.H[i] * wy.H[j] * wz.H[k];
                }
            }
        }
    }

    // Copy margins
    blitz::Range all = blitz::Range::all();
    blitz::Range left = blitz::Range(margin, 2 * margin - 1);
    blitz::Range left_margin = blitz::Range(0, margin - 1);
    blitz::Range right = blitz::Range(n_grid, n_grid + margin - 1);
    blitz::Range right_margin = blitz::Range(margin + n_grid, n_grid + 2 * margin - 1);
    // x-direction
    res_m(left, all, all) += res_m(right_margin, all, all);
    res_m(right, all, all) += res_m(left_margin, all, all);
    // y-direction
    res_m(all, left, all) += res_m(all, right_margin, all);
    res_m(all, right, all) += res_m(all, left_margin, all);
    // z-direction
    res_m(all, all, left) += res_m(all, all, right_margin);
    res_m(all, all, right) += res_m(all, all, left_margin);

    return res;
}

#endif // PROJECT_MASS_ASSIGNEMENT_H
