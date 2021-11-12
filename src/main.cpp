#include <chrono>

#include "aweights.h"
#include "blitz/array.h"
#include "io_utils.h"
#include "mass_assignement.h"
#include "tipsy.h"

int main() {
    using std::cout, std::endl, std::cerr;
    using namespace blitz;
    using namespace std::chrono;
    using real_t = float;

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

    auto n_grid = 64; // Number of grid cells per dimension.

    // 3rd order with margins
    auto start = chrono::high_resolution_clock::now();
    auto density_grid =
        assign_mass_with_margins<real_t, 3>(r, n_grid);

    // Turn density into overdensity.
    density_grid = density_grid / mean(density_grid) - 1;





    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "3rd order with margins: " << duration.count() << " milliseconds"
         << endl;
}
