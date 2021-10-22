#include <cstdint> 
#include "tipsy.h"
#include "aweights.hpp"
#include <chrono>

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

        float mass = 1.;
        
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

int write_to_csv(blitz::Array<float,2> array, const char* fname) {
    std::ofstream out;
    out.open(fname);
    if (!out.is_open()) {
        return 1;
    }

    for (auto i=0; i<array.extent(blitz::firstDim); ++i) {
        for (auto j=0; j<array.extent(blitz::secondDim); ++j) {
            out << array(i, j) << ",";
        }
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
    auto mass_grid = assign_mass<2>(r, n_grid, &wrap_modulo); 
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "2nd order wrap_modulo: " << duration.count() << " miliseconds" << endl;

    // 2nd order if-else
    start = chrono::high_resolution_clock::now();
    mass_grid = assign_mass<2>(r, n_grid, &wrap_if_else); 
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "2nd order if-else: " << duration.count() << " miliseconds" << endl;

    cout << endl;

    // 3nd order modulo
    start = chrono::high_resolution_clock::now();
    mass_grid = assign_mass<3>(r, n_grid, &wrap_modulo); 
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "3nd order wrap_modulo: " << duration.count() << " miliseconds" << endl;

    // 3nd order if-else
    start = chrono::high_resolution_clock::now();
    mass_grid = assign_mass<3>(r, n_grid, &wrap_if_else); 
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "3nd order if-else: " << duration.count() << " miliseconds" << endl;

    cout << endl;

    // 4nd order modulo
    start = chrono::high_resolution_clock::now();
    mass_grid = assign_mass<4>(r, n_grid, &wrap_modulo); 
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "4nd order wrap_modulo: " << duration.count() << " miliseconds" << endl;

    // 4nd order if-else
    start = chrono::high_resolution_clock::now();
    mass_grid = assign_mass<4>(r, n_grid, &wrap_if_else); 
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "4nd order if-else: " << duration.count() << " miliseconds" << endl;

}
