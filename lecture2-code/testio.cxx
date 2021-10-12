#include <cstdint> 
#include "tipsy.h"

/**
 * Return the nearest grid point for a particle .
 *
 * @param coords coordinates of the particle.
 * @param N size of the grid.
 */
blitz::TinyVector<int, 3>
ngp(blitz::TinyVector<float, 3> coords, int N) {
    blitz::TinyVector<int, 3> res(0);

    return blitz::ceil(N * (coords + 0.5) - 1);
}

/**
 * Return a grid of mass densities given a particle distribution.
 *
 * @param particles Particles to distribute
 * @param N Grid size in 3d.
 */
blitz::Array<float, 3> assign_mass_ngp(blitz::Array<float, 2> particles, int N) {
    blitz::Array<float, 3> res(N, N, N);
    for (auto particle=particles.begin(); particle!=particles.end(); ++particle) {
        float mass = 1.;
        auto grid_point = ngp(*particle, N);
        res(grid_point) += mass;
    }
    return res;
}

int main() {
    using std::cout, std::endl;
    using namespace blitz;

    // Load data.
    TipsyIO io;

    io.open("b0-final.std");
    int N = io.count();
    std::cout << "N = " << N << endl << endl;

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

    // Test ngp();
    TinyVector<float, 3> test_vector(0.);
    cout << "N/2: " << N/2 << endl;
    cout << "ngp of " << test_vector << ": " <<  ngp(test_vector, N) << endl;

    test_vector = -0.5, -0.499999, 0.5; 
    cout << "N: " << N << endl;
    cout << "ngp of " << test_vector << ": " <<  ngp(test_vector, N) << endl;
    
    cout << "element 0: " << r(0) << endl;
    cout << "ngp of element 0: " << ngp(r(0), N) << endl;

    cout << endl;

    // Test assign_mass_ngp()
    auto mass_grid = assign_mass_ngp(r, 400);
    int chunk = 5;
    cout << mass_grid(Range(0, chunk), Range(0, chunk), Range(0, chunk)) << endl;
}
