#include "blitz/array.h"
#include "mpi.h"
#include <complex>
#include <fftw3-mpi.h>
#include <fftw3.h>
#include <fstream>
#include <omp.h>
#include <stdio.h>
using namespace std;
using namespace blitz;

/** This examples transposes an N*N*N matrix along the first two axes with the
 * help of the fftw-mpi library. Run with mpi_size == 5 to get a nice example
 * with two slabs per  * rank, i.e. (2, N, N).
 *
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char *argv[]) {
    // MPI Init shizzle.
    int thread_support;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &thread_support);
    int mpi_size, mpi_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size); // Number of ranks
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank); // my rank
    fftw_mpi_init();

    const int N = 10;

    // Compute local sizes.
    ptrdiff_t local_n, local_start;
    ptrdiff_t local_n_transposed, local_start_transposed;
    ptrdiff_t alloc_local =
        fftw_mpi_local_size_3d_transposed(N, N, N, MPI_COMM_WORLD,
                               &local_n, &local_start,
                               &local_n_transposed, &local_start_transposed);

    assert(local_n == local_n_transposed);
    assert(local_start == local_start_transposed);

    // Uncomment to print local sizes.
    if (mpi_rank == 0) {
        printf("r, local_start, local_n, alloc_local: %d, %td, %td, %td\n", mpi_rank, local_start, local_n, alloc_local);
        printf("r, local_start, local_n, alloc_local: %d, %td, %td, %td\n", mpi_rank, local_start_transposed, local_n_transposed, alloc_local);
    }

    // Allocate local_nxNxN grid of doubles.
    Array<double, 3> grid(Range(local_start, local_start + local_n - 1),
                                   Range(0, N - 1),
                                   Range(0, N - 1));

    // Fill grid with numbers such that grid[x, y, z] = xyz where x, y, z are
    // the digits of the number, not a multiplication. This is only possible
    // because N=10 and allows to easily read of which axes where transposed.
    for (auto i = grid.begin(); i != grid.end(); ++i) {
        auto pos = i.position();
        *i = pos[0] * grid.stride(0) + pos[1] * grid.stride(1) + pos[2] * grid.stride(2);
    }

    // Print rank=0 grid.
    if (mpi_rank == 0) {
        cout << grid.stride() << endl;
        cout << "r=" << mpi_rank << ' ' << grid(Range::all(), Range::all(), Range::all()) << endl;
    }

    // Plan transpose.
    auto plan = fftw_mpi_plan_many_transpose(
        N, N, N, FFTW_MPI_DEFAULT_BLOCK,
        FFTW_MPI_DEFAULT_BLOCK, reinterpret_cast<double *>(grid.data()),
        reinterpret_cast<double *>(grid.data()), MPI_COMM_WORLD, FFTW_ESTIMATE);

    // Transpose and print.
    fftw_execute(plan);
    if (mpi_rank == 0)
        cout << grid(Range::all(), Range::all(), Range::all()) << endl;

    // Transpose back and print.
    fftw_execute(plan);
    if (mpi_rank == 0) {
        cout << grid(Range::all(), Range::all(), Range::all()) << endl;
    }

    // MPI cleanup shizzle.
    fftw_mpi_cleanup();
    MPI_Finalize();
}