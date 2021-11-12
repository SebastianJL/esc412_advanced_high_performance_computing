#include <chrono>

#include "aweights.h"
#include "blitz/array.h"
#include "mass_assignement.h"
#include "my_io.h"
#include "tipsy.h"

/*
 * Wrap the point i on a grid with size n_grid using the % operator.
 *
 * @param i Coordinate to wrap.
 * @param n_grid Wrap modulo n_grid.
 */
int wrap_modulo(int i, int n_grid) { return (i + n_grid) % n_grid; }

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
    using real_t = double;
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
    auto mass_grid_modulo = assign_mass<real_t, 2>(r, n_grid, &wrap_modulo);
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "2nd order wrap_modulo: " << duration.count() << " milliseconds"
         << endl;
    auto mass_grid_modulo_2d = project_3d_to_2d(mass_grid_modulo);
    write_to_csv(mass_grid_modulo_2d, "output/mass_grid_modulo_2d.csv");

    // 2nd order if-else
    start = chrono::high_resolution_clock::now();
    auto mass_grid_if_else = assign_mass<real_t, 2>(r, n_grid, &wrap_if_else);
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "2nd order if-else: " << duration.count() << " milliseconds"
         << endl;
    auto mass_grid_if_else_2d = project_3d_to_2d(mass_grid_if_else);
    write_to_csv(mass_grid_if_else_2d, "output/mass_grid_if_else_2d.csv");

    // 2nd order with margins
    start = chrono::high_resolution_clock::now();
    auto mass_grid_with_margins = assign_mass_with_margins<real_t, 2>(r, n_grid);
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
    mass_grid_modulo = assign_mass<real_t, 3>(r, n_grid, &wrap_modulo);
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "3nd order wrap_modulo: " << duration.count() << " milliseconds"
         << endl;

    // 3nd order if-else
    start = chrono::high_resolution_clock::now();
    mass_grid_if_else = assign_mass<real_t, 3>(r, n_grid, &wrap_if_else);
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "3nd order if-else: " << duration.count() << " milliseconds"
         << endl;

    // 3rd order with margins
    start = chrono::high_resolution_clock::now();
    mass_grid_with_margins = assign_mass_with_margins<real_t, 3>(r, n_grid);
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
    mass_grid_modulo = assign_mass<real_t, 4>(r, n_grid, &wrap_modulo);
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "4nd order wrap_modulo: " << duration.count() << " milliseconds"
         << endl;

    // 4nd order if-else
    start = chrono::high_resolution_clock::now();
    mass_grid_if_else = assign_mass<real_t, 4>(r, n_grid, &wrap_if_else);
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "4nd order if-else: " << duration.count() << " milliseconds"
         << endl;

    // 4th order with margins
    start = chrono::high_resolution_clock::now();
    mass_grid_with_margins = assign_mass_with_margins<real_t, 4>(r, n_grid);
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "4th order with margins: " << duration.count() << " milliseconds"
         << endl;
    cout << "modulo == if-else: " << all(mass_grid_modulo == mass_grid_if_else)
         << endl;
    cout << "modulo == margins: "
         << all(mass_grid_modulo == mass_grid_with_margins) << endl;
}
