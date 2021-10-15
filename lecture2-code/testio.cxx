#include <cstdint> 
#include "tipsy.h"

/**
 * Return the nearest grid point for a single coordinate.
 *
 * @param coord coordinate.
 * @param n_grid size of the grid.
 */
float ngp(float coord, int n_grid) {
    return std::ceil(n_grid * (coord + 0.5) - 1);
}

/**
 * Return a grid of mass densities given a particle distribution.
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

            // std::cout << *p << std::endl;
            // std::cout << p.position() << std::endl;
            // std::cout << particles(p.position(), blitz::Range::all()) << std::endl;
            // std::cout << grid_point << std::endl;
            // std::cout << std::endl;
        }
        float mass = 1.;
        res(grid_point) += mass;
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


int main() {
    using std::cout, std::endl;
    using namespace blitz;

    // Load data.
    TipsyIO io;

    io.open("b0-final.std");
    int N_grid = io.count();
    std::cout << "N_grid = " << N_grid << endl << endl;

    if (io.fail()) {
        std::cerr << "Unable to open file" << std::endl;
        abort();
    }

    Array<float,2> r(io.count(),3);
    io.load(r);
    
    // Print some data.
    cout << r(Range(0,10), Range::all()) << endl;

    cout << "max_x = " << max(r(Range::all(), 0)) << endl;
    cout << "max_y = " << max(r(Range::all(), 1)) << endl;
    cout << "max_z = " << max(r(Range::all(), 2)) << endl;

    cout << endl;

    cout << "min_x = " << min(r(Range::all(), 0)) << endl;
    cout << "min_y = " << min(r(Range::all(), 1)) << endl;
    cout << "min_z = " << min(r(Range::all(), 2)) << endl;

    cout << endl;

    // Test assign_mass_ngp()
    auto n_grid = 64; // Number of grid cells per dimension.
    auto mass_grid = assign_mass_ngp(r, n_grid); 
    int chunk = 5;
    cout << mass_grid(Range(0, chunk), Range(0, chunk), Range(0, chunk)) << endl;
    
    // Project to 2d grid.
    firstIndex i;
    secondIndex j;
    thirdIndex k;
    Array<float, 2> mass_grid_2d(n_grid, n_grid);
    mass_grid_2d = max(mass_grid, k);
    cout << mass_grid_2d << endl;

    // Write to file
    write_to_csv(mass_grid_2d, "mass_grid_2d.csv");
}
