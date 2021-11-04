#include <chrono>

#include "aweights.hpp"
#include "blitz/array.h"
#include "my_io.h"
#include "tipsy.h"

using real_t = double;

/**
 * Return the nearest grid point for a single coordinate.
 *
 * @param coord coordinate.
 * @param n_grid size of the grid.
 */
real_t ngp(real_t coord, int n_grid) {
    return std::floor(n_grid * (coord + 0.5));
}

/**
 * Return the grid coordinate from a euclidean coordinate.
 *
 * @param coord coordinate.
 * @param n_grid size of the grid.
 */
real_t grid_coordinate(real_t coord, int n_grid) {
    return n_grid * (coord + 0.5);
}

/**
 * Return a grid of mass densities given a particle distribution.
 *
 * The algorithm used is nearest grid point.
 *
 * @param particles Particles to distribute
 * @param n_grid Grid size in 3d.
 */
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

/**
 * Wrap the point i on a grid with size n_grid using the % operator.
 *
 * @param i Coordinate to wrap.
 * @param n_grid Wrap modulo n_grid.
 */
int wrap_modulo(int i, int n_grid) { return (i + n_grid) % n_grid; }

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
 * Return a grid of mass densities given a particle distribution.
 *
 * @param particles Particles to distribute
 * @param n_grid Grid size in 3d.
 * @param order Order of the mass assignment weights.
 */
template <int Order>
blitz::Array<real_t, 3> assign_mass(blitz::Array<real_t, 2> particles,
                                    int n_grid, int (*wrap)(int, int)) {
    blitz::Array<real_t, 3> res(n_grid, n_grid, n_grid);
    res = 0;
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
                    res(x, y, z) += wx.H[i] * wy.H[j] * wz.H[k];
                }
            }
        }
    }
    return res;
}

/**
 * Return a grid of mass densities given a particle distribution.
 *
 * Use a grid bigger with margins instead of using a modulo operation.
 *
 * @param particles Particles to distribute
 * @param n_grid Grid size in 3d.
 * @param order Order of the mass assignment weights.
 */
template <int Order>
blitz::Array<real_t, 3>
assign_mass_with_margins(blitz::Array<real_t, 2> particles, int n_grid) {
    using blitz::Array;
    using blitz::firstDim;
    using blitz::Range;

    int margin = (Order - 1);
    int n_grid_m = n_grid + 2 * margin;
    Array<real_t, 3> res_m(n_grid_m, n_grid_m, n_grid_m);
    res_m = 0;

    Range inside(Range(margin, n_grid + margin - 1));
    Array<real_t, 3> res = res_m(inside, inside, inside);

    for (auto row = 0; row < particles.extent(firstDim); ++row) {
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
    Range all = Range::all();
    Range left = Range(margin, 2 * margin - 1);
    Range left_margin = Range(0, margin - 1);
    Range right = Range(n_grid, n_grid + margin - 1);
    Range right_margin = Range(margin + n_grid, n_grid + 2 * margin - 1);
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

blitz::Array<real_t, 2>
project_3d_to_2d(const blitz::Array<real_t, 3> &grid_3d) {
    int n_grid = grid_3d.extent(blitz::firstDim);
    blitz::thirdIndex k;
    blitz::Array<real_t, 2> grid_2d(n_grid, n_grid);

    grid_2d = max(grid_3d, k);
    return grid_2d;
}

int main() {
    using std::cout, std::endl, std::cerr;
    using namespace blitz;
    using namespace std::chrono;

    // Load data.
    TipsyIO io;

    auto filename = "input/b0-final.std";
    io.open(filename);
    if (io.fail()) {
        cerr << "Unable to open file '" << filename << "'" << endl;
        abort();
    }

    uint64_t n_particle = io.count();
    cout << "n_particle = " << n_particle << endl;

    Array<real_t, 2> r(io.count(), 3);
    io.load(r);

    // Print some data.
    firstIndex i;
    secondIndex j;
    thirdIndex k;
    Array<real_t, 1> max_r(3);
    cout << "max_x = " << max(r(Range::all(), 0)) << endl;
    cout << "max_y = " << max(r(Range::all(), 1)) << endl;
    cout << "max_z = " << max(r(Range::all(), 2)) << endl;
    cout << "max = " << max(r(Range::all(), 2)) << endl;

    cout << endl;

    cout << "min_x = " << min(r(Range::all(), 0)) << endl;
    cout << "min_y = " << min(r(Range::all(), 1)) << endl;
    cout << "min_z = " << min(r(Range::all(), 2)) << endl;

    cout << endl;

    // Time different wrapping functions
    auto n_grid = 64; // Number of grid cells per dimension.

    // 2nd order modulo
    auto start = chrono::high_resolution_clock::now();
    auto mass_grid_modulo = assign_mass<2>(r, n_grid, &wrap_modulo);
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "2nd order wrap_modulo: " << duration.count() << " milliseconds"
         << endl;
    auto mass_grid_modulo_2d = project_3d_to_2d(mass_grid_modulo);
    write_to_csv(mass_grid_modulo_2d, "output/mass_grid_modulo_2d.csv");

    // 2nd order if-else
    start = chrono::high_resolution_clock::now();
    auto mass_grid_if_else = assign_mass<2>(r, n_grid, &wrap_if_else);
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "2nd order if-else: " << duration.count() << " milliseconds"
         << endl;
    auto mass_grid_if_else_2d = project_3d_to_2d(mass_grid_if_else);
    write_to_csv(mass_grid_if_else_2d, "output/mass_grid_if_else_2d.csv");

    // 2nd order with margins
    start = chrono::high_resolution_clock::now();
    auto mass_grid_with_margins = assign_mass_with_margins<2>(r, n_grid);
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "2nd order with margins: " << duration.count() << " milliseconds"
         << endl;
    auto mass_grid_with_margins_2d = project_3d_to_2d(mass_grid_with_margins);
    write_to_csv(mass_grid_with_margins_2d,
                 "output/mass_grid_with_margins_2d.csv");

    // Compare 2nd order methods.
    cout << "modulo == if-else: " << all(mass_grid_modulo == mass_grid_if_else)
         << endl;
    cout << "modulo == margins: "
         << all(mass_grid_modulo == mass_grid_with_margins) << endl;

    cout << endl;

    // 3nd order modulo
    start = chrono::high_resolution_clock::now();
    mass_grid_modulo = assign_mass<3>(r, n_grid, &wrap_modulo);
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "3nd order wrap_modulo: " << duration.count() << " milliseconds"
         << endl;

    // 3nd order if-else
    start = chrono::high_resolution_clock::now();
    mass_grid_if_else = assign_mass<3>(r, n_grid, &wrap_if_else);
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "3nd order if-else: " << duration.count() << " milliseconds"
         << endl;

    // 3rd order with margins
    start = chrono::high_resolution_clock::now();
    mass_grid_with_margins = assign_mass_with_margins<3>(r, n_grid);
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "3rd order with margins: " << duration.count() << " milliseconds"
         << endl;
    cout << "modulo == if-else: " << all(mass_grid_modulo == mass_grid_if_else)
         << endl;
    cout << "modulo == margins: "
         << all(mass_grid_modulo == mass_grid_with_margins) << endl;

    cout << endl;

    // 4nd order modulo
    start = chrono::high_resolution_clock::now();
    mass_grid_modulo = assign_mass<4>(r, n_grid, &wrap_modulo);
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "4nd order wrap_modulo: " << duration.count() << " milliseconds"
         << endl;

    // 4nd order if-else
    start = chrono::high_resolution_clock::now();
    mass_grid_if_else = assign_mass<4>(r, n_grid, &wrap_if_else);
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "4nd order if-else: " << duration.count() << " milliseconds"
         << endl;

    // 4th order with margins
    start = chrono::high_resolution_clock::now();
    mass_grid_with_margins = assign_mass_with_margins<4>(r, n_grid);
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "4th order with margins: " << duration.count() << " milliseconds"
         << endl;
    cout << "modulo == if-else: " << all(mass_grid_modulo == mass_grid_if_else)
         << endl;
    cout << "modulo == margins: "
         << all(mass_grid_modulo == mass_grid_with_margins) << endl;
}
