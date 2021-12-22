#include <fstream>
#include <complex>
#include <stdlib.h>
#include <sys/time.h>
#include <fftw3-mpi.h>
#include <omp.h>
#include <mpi.h>

#include "aweights.h"
#include "io_utils.h"
#include "tipsy.h"
#include "transformations.h"

using namespace std;

typedef double real_type;
typedef std::complex<real_type> complex_type;

typedef blitz::Array<real_type,2> array2D_r;
typedef blitz::Array<real_type,3> array3D_r;

typedef blitz::Array<complex_type,3> array3D_c;


// A simple timing function
double getTime() {
    struct timeval tv;
    struct timezone tz;
    gettimeofday(&tv, &tz);
    return tv.tv_sec + 1e-6*(double)tv.tv_usec;
}

// Read the particle file,
// return a 2d blitz array containing particle positions
array2D_r read_particles(string fname, int rank, int size){
    double t0, elapsed;
    TipsyIO io;

    io.open(fname);
    if (rank == 0) {
        cout << "(rank:" << rank << ") "
             << "Found " << io.count() << " particles." << endl;
    }
    if (io.fail()) {
        cerr << "Unable to open file" << endl;
        abort();
    }

    int N = io.count();                                  //Total particle number
    int n = (N+size-1) / size;                           //Particles per rank
    int iStart = rank*n;                                 //Start index
    int iEnd = iStart + n-1 < N ? iStart + n-1 : N-1;    //End index

    // Allocate the particle buffer
    array2D_r p(blitz::Range(iStart, iEnd), blitz::Range(0, 2));

    t0 = getTime();
    // Read the particles
    io.load(p);
    elapsed = getTime() - t0;

//    cout << "(rank:" << rank << ") " << "Particles read: " << p.length(0) << endl;
//    cout << "(rank:" << rank << ") " << "particle reading: " << elapsed << " s" << endl;
    return p;
}

// Write a blitz array in a csv file format
template<typename T>
void write_array(T A, const char* filename){
    cout << "Writing to " << filename << endl;
    ofstream ofs(filename);
    if (ofs.bad()){
        cerr << "Unable to write to file: " << filename << endl;
        abort();
    }

    ofs << A << endl;

    ofs.close();
}

// Projection of a 3D grid into a 2D grid (max pooling)
array2D_r project(const array3D_r& grid){
    auto shape = grid.shape();
    array2D_r ret(shape[0], shape[1]);
    blitz::thirdIndex k;
    ret = blitz::max(grid, k);
    return ret;
}


// Mass assignment for a single particle with order given by o
template<int Order>
void assign_mass(real_type x, real_type y, real_type z, array3D_r& grid){
    auto shape = grid.shape();
    auto N = shape[1];
    int i, j, k;
    AssignmentWeights<Order> wx((x + 0.5)*N);
    AssignmentWeights<Order> wy((y + 0.5)*N);
    AssignmentWeights<Order> wz((z + 0.5)*N);

    for (int ii = 0; ii < Order; ii++) {
        // In the x-direction the index has to be wrapped differently due to the margin.
        i = (wx.i+ii+N)%N;  // For max mpi_rank i in {grid.lbound(0)..N, 0 .. Order - 2}
        if (i < grid.lbound(0)) {
            i += N;
        }

        for (int jj = 0; jj < Order; jj++) {
            j = (wy.i + jj + N) % N;
            for (int kk = 0; kk < Order; kk++) {
                k = (wz.i + kk + N) % N;

                #pragma omp atomic
                grid(i, j, k) += wx.H[ii] * wy.H[jj] * wz.H[kk];
            }
        }
    }
}

// Mass assignment for a list of particles
template<int Order>
void assign_masses(array2D_r p, array3D_r &grid, int N_grid){
    int mpi_rank, mpi_size;

    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    double t0, elapsed;
    t0 = getTime();

    // Use a view of the grid without the fft padding
    array3D_r grid_nopad = grid(blitz::Range::all(),
                                blitz::Range::all(),
                                blitz::Range(0,N_grid - 1));
    #pragma omp parallel for  // NOLINT(openmp-use-default-none)
    for(auto i=p.lbound(0); i<=p.ubound(0); ++i){
        assign_mass<Order>(p(i, 0), p(i, 1), p(i, 2), grid_nopad);
    }

    // Exchange masses in margins.
    int dst = (mpi_rank + 1 + mpi_size) % mpi_size;
    int src = (mpi_rank - 1 + mpi_size) % mpi_size;

    if (dst != src){
        // Compute count and offset of the mass assignment margin in grid.
        auto shape = grid.shape();
        int count = 3*shape[1]*shape[2];
        int offset = (shape[0]-3)*shape[1]*shape[2];

        MPI_Request req;
        real_type *margin_recv = new real_type[count];
        MPI_Irecv(margin_recv, count, MPI_DOUBLE,
                  src, 0, MPI_COMM_WORLD, &req);

        MPI_Send(grid.dataFirst() + offset, count, MPI_DOUBLE,
                 dst, 0, MPI_COMM_WORLD);
        MPI_Wait(&req, MPI_STATUS_IGNORE);

        for(auto i=0; i<count; i++){
            grid.dataFirst()[i] += margin_recv[i]+1;
        }
        MPI_Barrier(MPI_COMM_WORLD); // Barrier just for time measurement
    }

    // Compute the average density per grid cell
    real_type local_mass = blitz::sum(grid_nopad);
    real_type total_mass = 0;

    MPI_Allreduce(&local_mass, &total_mass, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    real_type avg = total_mass / (N_grid*N_grid*N_grid);
    // Turn the density into the over-density.
    grid = (grid - avg) / avg;

    elapsed = getTime()-t0;
//    cout << "mass assignment: " << elapsed << " s" << endl;
}

void compute_fft(array3D_r grid, array3D_c fft_grid, int N, MPI_Comm comm){

    double t0, elapsed;

    // Create FFTW plan
    t0 = getTime();
    auto plan = fftw_mpi_plan_dft_r2c_3d(N,N,N,
            grid.dataFirst(),
            reinterpret_cast<fftw_complex*>(fft_grid.dataFirst()),
            comm,
            FFTW_ESTIMATE);
    elapsed = getTime()-t0;
//    cout << "fftw_plan creation: " << elapsed << " s" << endl;


    // Execute FFTW plan
    t0 = getTime();
    fftw_execute(plan);
    elapsed = getTime()-t0;
//    cout << "fftw_plan execution: " << elapsed << " s" << endl;

    // Destroy FFTW plan
    fftw_destroy_plan(plan);
}

void reduce_in_place(int rank, int root, void * buff, int size,
                     MPI_Datatype dtype, MPI_Op op, MPI_Comm comm){
    // In-place reduction with MPI (calls from rank 0 are different from others)
    if(rank==root){
        MPI_Reduce(MPI_IN_PLACE, buff, size,
                   dtype, op, root, comm);
    }
    else{
        MPI_Reduce(buff, NULL, size,
                   dtype, op, root, comm);
    }

}

void compute_pk(array3D_c fft_grid, int N){
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double t0, elapsed;
    int iNyquist = N / 2;
    int nBins = iNyquist;

    double k_max = sqrt(3) * (iNyquist+1);
    double scale = nBins * 1.0 / log(k_max);

    // LIN SPACED BINS:
    // bin_idx = floor( k_norm / k_max * nBins )

    // LOG SPACED BINS:
    // bin_idx = floor ( log(k_norm) * scale )
    //         = floor ( nBins * log(k_norm) / log(k_max) )

    blitz::Array<double,1>  log_k(nBins);
    blitz::Array<double,1>  power(nBins);
    blitz::Array<int64_t,1> n_power(nBins);


    // Fill arrays with zeros
    log_k = 0;
    power = 0;
    n_power = 0;

    // Mode ordering by fftw:
    // 0, 1, 2, ..., N/2, -(N/2-1), ..., -2, -1
    // 0, 1, 2, ..., N/2, -(N/2-1), ..., -2, -1
    // 0, 1, 2, ..., N/2

    t0 = getTime();
    for( auto index=fft_grid.begin(); index!=fft_grid.end(); ++index ) {
        auto pos = index.position();
        int kx = pos[0]>iNyquist ? N - pos[0] : pos[0];
        int ky = pos[1]>iNyquist ? N - pos[1] : pos[1];
        int kz = pos[2];

        if(kx==0 && ky==0 && kz==0) continue;

        int mult = (kz == 0) || (kz == iNyquist) ? 1 : 2;

        complex_type cplx_amplitude = *index;

        double k_norm = sqrt(kx*kx + ky*ky + kz*kz);

        int bin_idx = k_norm>0 ? floor(log(k_norm) * scale) : 0;

        assert(bin_idx>=0 && bin_idx<nBins);

        log_k(bin_idx)   += mult*log(k_norm);
        power(bin_idx)   += mult*norm(cplx_amplitude);
        n_power(bin_idx) += mult;
    }


    reduce_in_place(rank, 0, log_k.data(), nBins,
                     MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    reduce_in_place(rank, 0, power.data(), nBins,
                     MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    reduce_in_place(rank, 0, n_power.data(), nBins,
                     MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);

    elapsed = getTime() - t0;

//    cout << "P(k) measurement: " << elapsed << " s" << endl;

    if (rank == 0) {
        for (int i = 0; i < nBins; ++i) {
            if (n_power(i) > 0) {
                cout << exp(log_k(i) / n_power(i)) << " "
                     << power(i) / n_power(i) << endl;
            }
        }
    }
}


void finalize() {
    fftw_mpi_cleanup();
    MPI_Finalize();
}

/** Exchange particles such that each process only has the particles that fit
 * within its own part of the grid.
 *
 * This function assumes that each rank get's at least some part of the grid.
 * Meaning start != 0 iff rank != 0
 * @param particles
 * @param N_grid: Total size of the particle grid.
 * @param start: Local starting index in the grid. I.e. grid(i, 0, 0) is the
 * first element of the local grid.
 */
template<int Order, typename real_t>
void exchange_particles(blitz::Array<real_t, 2> &particles, int N_grid, ptrdiff_t start) {
    using blitz::Range;
    using blitz::Array;
    int mpi_rank, mpi_size;

    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    assert(mpi_size > 1);


    // Communicate domain decomposition with MPI_Allgather().
    blitz::Array<ptrdiff_t, 1> domain_boundaries(mpi_size + 1);
    domain_boundaries(mpi_size) = N_grid;
    MPI_Allgather(&start, 1, MPI_AINT, domain_boundaries.data(), 1, MPI_AINT, MPI_COMM_WORLD);

    // DEBUG. Check that the domain boundaries are strictly ascending.
    for (auto i=domain_boundaries.lbound(0) + 1; i<=domain_boundaries.ubound(0); ++i) {
        assert(domain_boundaries(i) > domain_boundaries(i-1));
    }


    // Count number of particles per slab.
    Array<int , 1> counts(N_grid); // Number of particles per slab.
    Array<int, 1> idx(Range(particles.lbound(0), particles.ubound(0))); // Slab index per particle.
    counts = 0;
    for (auto i = particles.lbound(0); i <= particles.ubound(0); ++i) {
        auto x_grid = grid_coordinate(particles(i, 0), N_grid);
        auto w = AssignmentWeights<Order>(x_grid);
        auto slab_index = wrap_if_else(w.i, N_grid);
        counts(slab_index) += 1;
        assert(idx.isInRange(i));
        idx(i) = slab_index;
    }

    // DEBUG Check that each particle has been assigned to exactly one slab.
    assert(sum(counts) == particles.length(0));
    for (auto i = counts.lbound(0); i <= counts.ubound(0); ++i) {
        assert(counts(i) >= 0);
    }

    // Sum particle counts.
    Array<int , 1> counts_t(N_grid); // Total counts up to slab per slab.
    // Instead of copying counts into counts_t it is referenced and then freed,
    // so counts "cannot be used afterwards".
    counts_t.reference(counts);
    counts.free();
    for (int i = counts_t.lbound(0) + 1; i <= counts_t.ubound(0); ++i) {
        counts_t(i) += counts_t(i-1);
    }

    // DEBUG Check that counts_t have been correctly summed and are ascending.
    assert(counts_t(counts_t.ubound(0)) == particles.length(0));
    for (auto i = counts_t.lbound(0); i <= counts_t.ubound(0); ++i) {
        assert(counts_t(i) >= 0);
    }


    // Infer counts/offsets per rank from counts_t.
    blitz::Array<int , 1> offsets_per_rank(mpi_size); // Offset of each rank in particles array.
    blitz::Array<int , 1> counts_per_rank(mpi_size); // Number of particles per rank.
    counts_per_rank = 0;
    counts_per_rank(0) = counts_t(domain_boundaries(1) - 1);
    for (auto i = 1; i < mpi_size; ++i) {
        assert(counts_per_rank.isInRange(i));
        counts_per_rank(i) = counts_t(domain_boundaries(i+1) - 1) - counts_t(domain_boundaries(i) - 1);
    }
    offsets_per_rank = 0;
    for (auto i = 1; i < mpi_size; ++i) {
        offsets_per_rank(i) = offsets_per_rank(i-1) + counts_per_rank(i-1);
    }

    // DEBUG.
    assert(offsets_per_rank(offsets_per_rank.lbound(0)) == 0);
    assert(offsets_per_rank(offsets_per_rank.ubound(0)) <= particles.length(0));
    assert(sum(counts_per_rank) == particles.length(0));
    for (auto i = counts_per_rank.lbound(0); i <= counts_per_rank.ubound(0); ++i) {
        assert(counts_per_rank(i) >= 0);
    }


    // Order particles by slab (out-of-place).
    // p_out does only have the extent() of particles but not the lbound().
    // The reason is, that counts_t only counts the offset from the start of the
    // particle array and doesn't consider the lbound().
    Array<real_t, 2> p_out(particles.extent());
    auto all = blitz::Range::all();

    for (int i = particles.ubound(0); i >= particles.lbound(0); --i){
        assert(p_out.isInRange(counts_t(idx(i)) - 1));
        p_out(counts_t(idx(i)) - 1, all) = particles(i, all);
        counts_t(idx(i)) -= 1;
    }
    // Set p_out.lbound() to particles.lbound() such that particles keeps it's
    // lbound() when it references p_out.
    p_out.reindexSelf(particles.lbound());
    particles.reference(p_out);
    p_out.free();


    // DEBUG Check that particles are in order.
    auto current_index = 0;
    for (auto i = particles.lbound(0); i <= particles.ubound(0); ++i) {
        auto x_grid = grid_coordinate(particles(i, 0), N_grid);
        auto w = AssignmentWeights<Order>(x_grid);
        auto slab_index = wrap_if_else(w.i, N_grid);
        assert(slab_index >= current_index);
        current_index = slab_index;
    }


    // Communicate counts_per_rank with MPI_Alltoall()
    Array<int, 1> counts_from_rank(mpi_size);
    counts_from_rank = 0;
    MPI_Alltoall(counts_per_rank.data(), 1, MPI_INT,
                 counts_from_rank.data(), 1, MPI_INT, MPI_COMM_WORLD);
    Array<int, 1> offsets_from_rank(mpi_size);
    offsets_from_rank(0) = 0;
    for (auto i = counts_from_rank.lbound(0) + 1; i <= counts_from_rank.ubound(0); ++i) {
        offsets_from_rank(i) = counts_from_rank(i-1) + offsets_from_rank(i-1);
    }


    // Exchange particles with MPI_Alltoallv().
    Array<double, 2> particles_balanced(sum(counts_from_rank), 3);
    particles_balanced = 0;

    MPI_Datatype dt_particle;
    MPI_Type_contiguous(3, MPI_DOUBLE, &dt_particle);
    MPI_Type_commit(&dt_particle);

    MPI_Alltoallv(particles.data(), counts_per_rank.data(),
                  offsets_per_rank.data(), dt_particle,
                  particles_balanced.data(), counts_from_rank.data(),
                  offsets_from_rank.data(), dt_particle,
                  MPI_COMM_WORLD);

    MPI_Type_free(&dt_particle);
    particles.reference(particles_balanced);
    particles_balanced.free();

    // DEBUG check that all particles now belong into their own rank.
    for (auto i = particles.lbound(0); i <= particles.ubound(0); ++i) {
        auto x_grid = grid_coordinate(particles(i, 0), N_grid);
        auto w = AssignmentWeights<Order>(x_grid);
        auto slab_index = wrap_modulo(w.i, N_grid);
        assert(slab_index >= domain_boundaries(mpi_rank));
        assert(slab_index < domain_boundaries(mpi_rank+1));
    }
}

int main(int argc, char *argv[]) {
    using namespace blitz;
    static const int Order = 4;
    int thread_support;
    int mpi_rank, mpi_size;

    MPI_Init_thread(&argc, &argv,
        MPI_THREAD_FUNNELED,
        &thread_support);

    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    fftw_mpi_init();

    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());

    int N_grid = 64;
    string fname = "input/b0-final.std";
    array2D_r p = read_particles(fname, mpi_rank, mpi_size);

    // Compute local sizes
    ptrdiff_t alloc_local, local_n0, local_0_start;
    alloc_local = fftw_mpi_local_size_3d(N_grid, N_grid, N_grid/2 + 1, MPI_COMM_WORLD,
                                         &local_n0, &local_0_start);
    assert(local_n0 >= Order - 1); // Check that particles cannot be spread over more then two ranks.

    // Exchange particles
    if (mpi_size > 1) {
        if (local_n0 == 0) {
            cout << "(rank: " << mpi_rank << ") local_n0: " << local_n0 << endl;
        }
        assert(local_n0 != 0);
        exchange_particles<Order>(p, N_grid, local_0_start);
    } else {
        assert(mpi_size == 1);
    }

    // Allocate overlapping grids for mass assignment and fft.
    alloc_local += (Order - 1) * N_grid * (N_grid/2 + 1);  // Add margin for overhanging particles to alloc local.
    auto buffer = fftw_alloc_real(2 * alloc_local);
    GeneralArrayStorage<3> storage;
    storage.base() = local_0_start, 0, 0;
    // Grid containing mass assignment margin and fft padding.
    array3D_r raw_grid(
        buffer,
        shape(local_n0 + Order - 1, N_grid, 2 * (N_grid / 2 + 1)),
        neverDeleteData,
        storage);
    // Grid for real fft input.
    array3D_r grid_r(
        buffer,
        shape(local_n0, N_grid, 2 * (N_grid / 2 + 1)),
        neverDeleteData,
        storage);
    // Grid for complex fft output.
    array3D_c grid_c(
        reinterpret_cast<complex_type*> (buffer),
        shape(local_n0, N_grid, (N_grid / 2 + 1)),
        neverDeleteData,
        storage);


    assign_masses<Order>(p, raw_grid, N_grid);
    p.free();

    // Compute the fft of the over-density field
    compute_fft(grid_r, grid_c, N_grid, MPI_COMM_WORLD);

    // Compute the power spectrum
    compute_pk(grid_c, N_grid);

    finalize();
    return 0;
}
