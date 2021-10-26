#include <cstdint> 
#include <chrono>
#include "blitz/array.h"
#include "tipsy.h"
#include "aweights.hpp"


/**
 * Return the nearest grid point for a single coordinate.
 *
 * @param coord coordinate.
 * @param n_grid size of the grid.
 */
float ngp(float coord, int n_grid) {
    return std::floor(n_grid * (coord + 0.5));
}


/**
 * Return the grid coordinate from a euclidian coordinate. 
 *
 * @param coord coordinate.
 * @param n_grid size of the grid.
 */
float grid_coordinate(float coord, int n_grid) {
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
blitz::Array<float, 3>
assign_mass_ngp(blitz::Array<float, 2> particles, int n_grid) {
    blitz::Array<float, 3> res(n_grid, n_grid, n_grid);
    for (auto i=0; i<particles.extent(blitz::firstDim); ++i) {
        blitz::TinyVector<float, 3> grid_point;
        for (auto j=0; j<particles.extent(blitz::secondDim); ++j) {
            float p = particles(i, j);
            grid_point(j) = ngp(p, n_grid);
        }
        float mass = 1.;
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
int wrap_modulo(int i, int n_grid) {
    return (i + n_grid) % n_grid;
}

/**
 * Wrap the point i on a grid with size n_grid using if else statements.
 *
 * @param i Coordinate to wrap.
 * @param n_grid Wrap modulo n_grid.
 */
int wrap_if_else(int i, int n_grid) {
    if (i >= n_grid) return i - n_grid;
    if (i < 0) return i + n_grid;
    return i;
}

/**
 * Return a grid of mass densities given a particle distribution.
 *
 * @param particles Particles to distribute
 * @param n_grid Grid size in 3d.
 * @param order Order of the mass asignement weights.
 */
template<int Order>
blitz::Array<float, 3>
assign_mass(blitz::Array<float, 2> particles, int n_grid, int (*wrap) (int, int)) {
    blitz::Array<float, 3> res(n_grid, n_grid, n_grid);
    res = 0;
    for (auto row=0; row<particles.extent(blitz::firstDim); ++row) {
        float p = particles(row, 0);
        float p_grid = grid_coordinate(p, n_grid);
        AssignmentWeights<Order> wx(p_grid);

        p = particles(row, 1);
        p_grid = grid_coordinate(p, n_grid);
        AssignmentWeights<Order> wy(p_grid);

        p = particles(row, 2);
        p_grid = grid_coordinate(p, n_grid);
        AssignmentWeights<Order> wz(p_grid);

        for (auto i=0; i<Order; ++i) {
            for (auto j=0; j<Order; ++j) {
                for (auto k=0; k<Order; ++k) {
                    int x = wrap(wx.i + i, n_grid);
                    int y = wrap(wy.i + i, n_grid);
                    int z = wrap(wz.i + i, n_grid);
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
 * @param order Order of the mass asignement weights.
 */
template<int Order>
blitz::Array<float, 3>
assign_mass_with_margins(blitz::Array<float, 2> particles, int n_grid) {
    using blitz::Array;
    using blitz::Range;
    using blitz::firstDim;
    
    int margin = (Order - 1);
    Array<float, 3> res_with_margins(
            Range(0, n_grid + 2 * margin - 1),
            Range(0, n_grid + 2 * margin - 1),
            Range(0, n_grid + 2 * margin - 1)
            );
    Array<float, 3> res = res_with_margins(
            Range(margin, n_grid + margin - 1),
            Range(margin, n_grid + margin - 1),
            Range(margin, n_grid + margin - 1)
            );
    res = 0;
    for (auto row=0; row<particles.extent(firstDim); ++row) {
        float p = particles(row, 0);
        float p_grid = grid_coordinate(p, n_grid);
        AssignmentWeights<Order> wx(p_grid);

        p = particles(row, 1);
        p_grid = grid_coordinate(p, n_grid);
        AssignmentWeights<Order> wy(p_grid);

        p = particles(row, 2);
        p_grid = grid_coordinate(p, n_grid);
        AssignmentWeights<Order> wz(p_grid);

        for (auto i=0; i<Order; ++i) {
            for (auto j=0; j<Order; ++j) {
                for (auto k=0; k<Order; ++k) {
                    int x = wx.i + i;
                    int y = wy.i + i;
                    int z = wz.i + i;
                    res(x, y, z) += wx.H[i] * wy.H[j] * wz.H[k];
                }
            }
        }

    }
    // Copy margins
    // x-direction
    Range all = Range::all();
    res(Range(0, margin), all, all) += res(Range(n_grid, n_grid + margin - 1), all, all);
    res(Range(n_grid - margin, n_grid - 1), all, all) += res(Range(-margin, - 1), all, all);
    // y-direction
    res(all, Range(0, margin), all) += res(all, Range(n_grid, n_grid + margin - 1), all);
    res(all, Range(n_grid - margin, n_grid - 1), all) += res(all, Range(-margin, - 1), all);
    // z-direction
    res(all, all, Range(0, margin)) += res(all, all, Range(n_grid, n_grid + margin - 1));
    res(all, all, Range(n_grid - margin, n_grid - 1)) += res(all, all, Range(-margin, - 1));
    return res;
}

int write_to_csv(blitz::Array<float,2> array, const char* fname) {
    std::ofstream out;
    out.open(fname);
    if (!out.is_open()) {
        return 1;
    }
    
    int i_max = array.extent(blitz::firstDim);
    int j_max = array.extent(blitz::secondDim);
    for (auto i=0; i<i_max; ++i) {
        for (auto j=0; j<j_max - 1; ++j) {
            out << array(i, j) << ",";
        }
        out << array(i, j_max -1);
        out << std::endl;
    }
    out.close();
    return 0;
}

blitz::Array<float, 2> project_3d_to_2d(blitz::Array<float, 3> grid_3d) {
    int n_grid = grid_3d.extent(blitz::firstDim);
    blitz::thirdIndex k;
    blitz::Array<float, 2> grid_2d(n_grid, n_grid);

    grid_2d =  max(grid_3d, k);
    return grid_2d;
}


int main() {
    using std::cout, std::endl;
    using namespace blitz;
    using namespace std::chrono;

    // Load data.
    TipsyIO io;

    io.open("b0-final.std");
    int n_particle = io.count();
    std::cout << "n_particle = " << n_particle << endl << endl;

    if (io.fail()) {
        std::cerr << "Unable to open file" << std::endl;
        abort();
    }

    Array<float,2> r(io.count(),3);
    io.load(r);

    // Print some data.
    firstIndex i;
    secondIndex j;
    thirdIndex k;
    Array<float, 1> max_r(3);
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
    cout << "2nd order wrap_modulo: " << duration.count() << " miliseconds" << endl;
    auto mass_grid_modulo_2d = project_3d_to_2d(mass_grid_modulo);

    // Write to file
    write_to_csv(mass_grid_modulo_2d, "mass_grid_modulo_2d.csv");

    // 2nd order if-else
    start = chrono::high_resolution_clock::now();
    auto mass_grid_if_else = assign_mass<2>(r, n_grid, &wrap_if_else); 
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "2nd order if-else: " << duration.count() << " miliseconds" << endl;

    // 2nd order with margins
    start = chrono::high_resolution_clock::now();
    auto mass_grid_with_margins = assign_mass_with_margins<2>(r, n_grid); 
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "2nd order with margins: " << duration.count() << " miliseconds" << endl;
    cout << "modulo == if-else: " << all(mass_grid_modulo == mass_grid_if_else) << endl;
    cout << "modulo == margins: " << all(mass_grid_modulo == mass_grid_with_margins) << endl;
    auto mass_grid_with_margins_2d = project_3d_to_2d(mass_grid_with_margins);

    // Write to file
    write_to_csv(mass_grid_with_margins_2d, "mass_grid_with_margins_2d.csv");

    cout << endl;

    // 3nd order modulo
    start = chrono::high_resolution_clock::now();
    mass_grid_modulo = assign_mass<3>(r, n_grid, &wrap_modulo); 
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "3nd order wrap_modulo: " << duration.count() << " miliseconds" << endl;

    // 3nd order if-else
    start = chrono::high_resolution_clock::now();
    mass_grid_if_else = assign_mass<3>(r, n_grid, &wrap_if_else); 
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "3nd order if-else: " << duration.count() << " miliseconds" << endl;

    // 3rd order with margins
    start = chrono::high_resolution_clock::now();
    mass_grid_with_margins = assign_mass_with_margins<3>(r, n_grid); 
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "3rd order with margins: " << duration.count() << " miliseconds" << endl;
    cout << "modulo == if-else: " << all(mass_grid_modulo == mass_grid_if_else) << endl;
    cout << "modulo == margins: " << all(mass_grid_modulo == mass_grid_with_margins) << endl;

    cout << endl;

    // 4nd order modulo
    start = chrono::high_resolution_clock::now();
    mass_grid_modulo = assign_mass<4>(r, n_grid, &wrap_modulo); 
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "4nd order wrap_modulo: " << duration.count() << " miliseconds" << endl;

    // 4nd order if-else
    start = chrono::high_resolution_clock::now();
    mass_grid_if_else = assign_mass<4>(r, n_grid, &wrap_if_else); 
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "4nd order if-else: " << duration.count() << " miliseconds" << endl;

    // 4th order with margins
    start = chrono::high_resolution_clock::now();
    mass_grid_with_margins = assign_mass_with_margins<4>(r, n_grid); 
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "4th order with margins: " << duration.count() << " miliseconds" << endl;
    cout << "modulo == if-else: " << all(mass_grid_modulo == mass_grid_if_else) << endl;
    cout << "modulo == margins: " << all(mass_grid_modulo == mass_grid_with_margins) << endl;

}
