
#include "transformations.h"
#include "io_utils.h"

using namespace std;

int main() {

    // test grid_coordinate()
    int n_grid = 64;
    float x_coordinate = 0;
    float g_coord = grid_coordinate(x_coordinate, n_grid);
    cout << "grid_coordinate(" << x_coordinate << ", " << n_grid
         << ") = " << g_coord << endl;
    assert(g_coord == 32);

    n_grid = 64;
    x_coordinate = -0.5;
    g_coord = grid_coordinate(x_coordinate, n_grid);
    cout << "grid_coordinate(" << x_coordinate << ", " << n_grid
         << ") = " << g_coord << endl;
    assert(g_coord == 0);

    n_grid = 64;
    x_coordinate = 0.5;
    g_coord = grid_coordinate(x_coordinate, n_grid);
    cout << "grid_coordinate(" << x_coordinate << ", " << n_grid
         << ") = " << g_coord << endl;
    assert(g_coord == 64);

    // test particle_rank()
    const int Order = 4;
    g_coord = 30;
    const int len = 4;
    ptrdiff_t starting_indices[len] = {0, 16, 32, 48};
    int p_rank = particle_rank<Order>(g_coord, starting_indices, len);
    cout << "particle_rank<" << Order << ">(" << g_coord << ", " << n_grid << ", "
         << sprint_array(starting_indices, len) << ", "
         << len << ") = " << p_rank << endl;
    assert(p_rank == 1);

    float g_coords[] = {1, 1.4, 1.5, 2};
    for (float g_coord : g_coords) {
        p_rank = particle_rank<4>(g_coord, starting_indices, len);
        cout << "particle_rank<" << Order << ">(" << g_coord << ", " << n_grid << ", "
             << sprint_array(starting_indices, len) << ", " << len
             << ") = " << p_rank << endl;
        if (g_coord < 1.5) {
            assert(p_rank == 3);
        }
        else {
            assert(p_rank == 0);
        }
    }


}