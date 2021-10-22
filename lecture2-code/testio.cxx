#include <cstdint> 
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
 * Return a grid of mass densities given a particle distribution.
 *
 * @param particles Particles to distribute
 * @param n_grid Grid size in 3d.
 * @param order Order of the mass asignement weights.
 */
template<int Order>
blitz::Array<float, 3>
assign_mass(blitz::Array<float, 2> particles, int n_grid) {
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
                    int x = (wx.i + i + n_grid) % n_grid;
                    int y = (wy.i + j + n_grid) % n_grid;
                    int z = (wz.i + k + n_grid) % n_grid;
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
    auto mass_grid_ngp = assign_mass_ngp(r, n_grid); 
    auto mass_grid_2nd = assign_mass<2>(r, n_grid); 
    auto mass_grid_3rd = assign_mass<3>(r, n_grid); 
    auto mass_grid_4th = assign_mass<4>(r, n_grid); 
    
    // Project to 2d grid.
    auto mass_grid_ngp_2d = project_3d_to_2d(mass_grid_ngp);
    auto mass_grid_2nd_2d = project_3d_to_2d(mass_grid_2nd);
    auto mass_grid_3rd_2d = project_3d_to_2d(mass_grid_3rd);
    auto mass_grid_4th_2d = project_3d_to_2d(mass_grid_4th);

    // Write to file
    write_to_csv(mass_grid_ngp_2d, "mass_grid_ngp_2d.csv");
    write_to_csv(mass_grid_2nd_2d, "mass_grid_2nd_2d.csv");
    write_to_csv(mass_grid_3rd_2d, "mass_grid_3rd_2d.csv");
    write_to_csv(mass_grid_4th_2d, "mass_grid_4th_2d.csv");


}
