#include <chrono>
#include <fftw3.h>

#include "aweights.h"
#include "blitz/array.h"
#include "io_utils.h"
#include "mass_assignement.h"
#include "tipsy.h"

/*
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

int main() {
    using std::cout, std::endl, std::cerr;
    using namespace blitz;
    using namespace std::chrono;
    using real_t = float;

    // Load data.
    TipsyIO io;
    auto infile = "input/b0-final.std";
    io.open(infile);
    if (io.fail()) {
        cerr << "Unable to open file '" << infile << "'" << endl;
        abort();
    }
    uint64_t n_particle = io.count();
    cout << "n_particle = " << n_particle << endl;
    Array<real_t, 2> r(io.count(), 3);
    io.load(r);

    // Start clock
    auto start = chrono::high_resolution_clock::now();

    // Compute fft
    auto n_grid = 64; // Number of grid cells per dimension.
    auto *buffer =
        (real_t *)fftw_malloc(n_grid * n_grid * n_grid * sizeof(real_t));
    Array<real_t, 3> grid(buffer, shape(n_grid, n_grid, n_grid),
                          neverDeleteData);
    assign_mass<real_t, 3>(r, n_grid, wrap_if_else, grid);

    grid = grid / mean(grid) - 1;  // Overdensity

    // Stop clock
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "3rd order if-else: " << duration.count() << " milliseconds"
         << endl;

    // Output result
    auto grid_2d = project_3d_to_2d(grid);
    auto outfile = "output/fourier_transform.csv";
    write_to_csv(grid_2d, outfile);
}
