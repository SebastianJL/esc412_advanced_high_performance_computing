#include <fstream>
#include <typeinfo>
#include <complex>
#include <stdlib.h>
#include <sys/time.h>
#include <fftw3-mpi.h>
#include <omp.h>
#include <mpi.h>

#include "aweights.hpp"
#include "tipsy.h"

using namespace std;

typedef double real_type;
typedef std::complex<real_type> complex_type;

typedef blitz::Array<real_type,2> array2D_r;
typedef blitz::Array<real_type,3> array3D_r;

typedef blitz::Array<complex_type,3> array3D_c;

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

    if(rank==0)
        cout << "Found "<<io.count() << " particles."<<endl;

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

    MPI_Barrier(MPI_COMM_WORLD); // Barrier just for time measurement
    elapsed = getTime() - t0;
    if(rank==0){
        cout << "particle reading: " << elapsed << " s" << endl;
    }
    return p;
}


// Write a blitz array in a csv file format
template<typename T>
void write_array(T A, const char* filename){
    cout << "Writing to " << filename << endl;
    ofstream ofs(filename);
    if (ofs.bad()){
        cerr << "Unable to write to file: " << filename << endl;
        abort();;
    }

    ofs << A << endl;

    ofs.close();
    return;
}

// Projection of a 3D grid into a 2D grid (max pooling)
array2D_r project(array3D_r grid){
    auto shape = grid.shape();
    array2D_r ret(shape[0], shape[1]);
    blitz::thirdIndex k;
    ret = blitz::max(grid, k);
    return ret;
}


// Mass assignment for a single particle with order given by o
template<int o>
void _assign_mass(real_type x, real_type y, real_type z, array3D_r& grid){
    auto shape = grid.shape();
    auto N = shape[1];
    int i, j, k;
    AssignmentWeights<o> wx((x + 0.5)*N);
    AssignmentWeights<o> wy((y + 0.5)*N);
    AssignmentWeights<o> wz((z + 0.5)*N);

    for(int ii=0; ii<o; ii++){
        i = (wx.i+ii+N)%N;

        if(i<grid.lbound(0)){
            i = grid.ubound(0)-2 + i;
        }

        for(int jj=0; jj<o; jj++){
            j = (wy.i+jj+N)%N;
            for(int kk=0; kk<o; kk++){
                k = (wz.i+kk+N)%N;

                #pragma omp atomic
                grid(i,j,k) += wx.H[ii]*wy.H[jj]*wz.H[kk];
            }
        }
    }
    return;
}

// Wrapper for templated mass assignment
void assign_mass(int o, real_type x, real_type y, real_type z, array3D_r& grid){
    switch(o){
        case 1: _assign_mass<1>(x,y,z,grid); break;
        case 2: _assign_mass<2>(x,y,z,grid); break;
        case 3: _assign_mass<3>(x,y,z,grid); break;
        case 4: _assign_mass<4>(x,y,z,grid); break;
        default:
            cerr << "Incorrect mass assignment order: " << o << endl;
            abort();
    }
}

// Mass assignment for a list of particles
void assign_masses(int o, array2D_r p, array3D_r &grid, int rank, int size){
    double t0, elapsed;
    t0 = getTime();

    auto shape = grid.shape();
    auto N = shape[1];

    // Use a view of the grid without the padding
    array3D_r grid_nopad = grid(blitz::Range::all(),
                                blitz::Range::all(),
                                blitz::Range(0,N-1));

    #pragma omp parallel for
    for(auto i=p.lbound(0); i<=p.ubound(0); ++i){
        assign_mass(o, p(i,0), p(i,1), p(i,2), grid_nopad);
    }
    real_type avg = blitz::sum(grid) / (N*N*N);

    MPI_Allreduce(MPI_IN_PLACE, &avg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // Turn the density into the over-density
    grid = (grid - avg) / avg;

    int dst = (rank+1+size)%size;
    int src = (rank-1+size)%size;

    if (dst != src){
        int count = 3*shape[1]*shape[2];
        int offset = (shape[0]-3)*shape[1]*shape[2];

        MPI_Request req;
        real_type *data = new real_type[count];
        MPI_Irecv(data, count, MPI_DOUBLE,
                  src, 0, MPI_COMM_WORLD, &req);

        MPI_Send(grid.dataFirst() + offset, count, MPI_DOUBLE,
                  dst, 0, MPI_COMM_WORLD);
        MPI_Wait(&req, MPI_STATUS_IGNORE);

        for(auto i=0; i<count; i++){
            grid.dataFirst()[i] += data[i]+1;
        }
        MPI_Barrier(MPI_COMM_WORLD); // Barrier just for time measurement
    }

    elapsed = getTime()-t0;
    if(rank==0){
        cout << "mass assignment: " << elapsed << " s" << endl;
    }
}

void compute_fft(array3D_r grid, array3D_c fft_grid, int N, MPI_Comm comm){
    int rank;
    MPI_Comm_rank(comm, &rank);
    double t0, elapsed;

    // Create FFTW plan
    t0 = getTime();
    auto plan = fftw_mpi_plan_dft_r2c_3d(N,N,N,
            grid.dataFirst(),
            reinterpret_cast<fftw_complex*>(fft_grid.dataFirst()),
            comm,
            FFTW_ESTIMATE);

    MPI_Barrier(MPI_COMM_WORLD); // Barrier just for time measurement
    elapsed = getTime()-t0;

    if(rank==0){
        cout << "fftw_plan creation: " << elapsed << " s" << endl;
    }


    // Execute FFTW plan
    t0 = getTime();
    fftw_execute(plan);

    elapsed = getTime()-t0;
    if(rank==0){
        cout << "fftw_plan execution: " << elapsed << " s" << endl;
    }

    // Destroy FFTW plan
    fftw_destroy_plan(plan);
}

void compute_pk(array3D_c fft_grid, int N, int rank){
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
    for(auto index=fft_grid.begin(); index!=fft_grid.end(); ++index){
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
    if(rank==0){
        cout << "P(k) measurement: " << elapsed << " s" << endl;

        for(int i=0; i<nBins; ++i) {
            if (n_power(i)>0) {
                cout << exp(log_k(i)/n_power(i)) << " " << power(i)/n_power(i) <<" "<< n_power(i)<< endl;
            }
        }
    }

}

int slab_idx(real_type pos_x, int N, int* domain, int size){
    if(size==0) return 0;
    AssignmentWeights<4> wx((pos_x+0.5)*N);
    int i = (wx.i+N)%N;
    int idx;
    for(idx=0; idx<size-1; idx++){
        if(i < domain[idx+1])
            break;
    }
    return idx;
}

void sort_particles(array2D_r &p, int N, int* domain, int size){
    int iStart = p.lbound(0);
    int iEnd   = p.ubound(0);
    vector<int> count_list(size+1, 0);
    array2D_r out(blitz::Range(iStart, iEnd), blitz::Range(0, 2));

    int idx;
    #pragma omp parallel for
    for (int i=iStart; i<=iEnd; i++){
        idx = slab_idx(p(i,0), N, domain, size);
        #pragma omp atomic
        count_list[idx] ++;
    }

    for (int i=1; i<count_list.size(); i++)
        count_list[i] += count_list[i-1];

    int pos;

    #pragma omp parallel for
    for(int i=iEnd;i>= iStart; i--){
        idx = slab_idx(p(i,0), N, domain, size);

        #pragma omp atomic capture
        {pos = count_list[idx]; count_list[idx]--;}

        out(iStart + pos-1, 0) = p(i, 0);
        out(iStart + pos-1, 1) = p(i, 1);
        out(iStart + pos-1, 2) = p(i, 2);
    }

    #pragma omp parallel for
    for (int i=iStart; i<=iEnd; i++){
        p(i,0) = out(i,0);
        p(i,1) = out(i,1);
        p(i,2) = out(i,2);
    }
}


array2D_r balance_particles(array2D_r p, ptrdiff_t local_0_start,
                            int N, int rank, int size){
    double t0, elapsed;
    t0 = getTime();
    int* domain = new int[size];

    MPI_Allgather(&local_0_start, 1, MPI_INTEGER,
            domain, 1, MPI_INTEGER, MPI_COMM_WORLD);

    sort_particles(p, N, domain, size);

    int* send_counts = new int[size];
    int* recv_counts = new int[size];
    int* send_offset = new int[size];
    int* recv_offset = new int[size];

    for(auto i=0; i<size; i++) {
        send_counts[i]=0;
        send_offset[i]=0;
        recv_offset[i]=0;
    }

    #pragma omp parallel for
    for(auto i=p.lbound(0); i<=p.ubound(0); ++i){
        int idx = slab_idx(p(i,0), N, domain, size);
        #pragma omp atomic
        send_counts[idx]++;
    }

    MPI_Alltoall(send_counts, 1, MPI_INTEGER, recv_counts,1, MPI_INTEGER, MPI_COMM_WORLD);

    int n_part = 0;

    for(int i=0; i<size; i++){
        n_part += recv_counts[i];
        recv_counts[i] *= 3;
        send_counts[i] *= 3;
        if(i>=1){
            send_offset[i] = send_offset[i-1] + send_counts[i-1];
            recv_offset[i] = recv_offset[i-1] + recv_counts[i-1];
        }
    }

    array2D_r p_balanced(n_part, 3);

    MPI_Alltoallv(p.dataFirst(), send_counts,
                  send_offset, MPI_DOUBLE, p_balanced.dataFirst(),
                  recv_counts, recv_offset, MPI_DOUBLE,
                  MPI_COMM_WORLD);

    elapsed = getTime() - t0;
    if(rank==0)
        cout << "load balancing: " << elapsed << " s" << endl;

    return p_balanced;
}



int main(int argc, char *argv[]) {
    int thread_support;
    int rank, size;

    MPI_Init_thread(&argc, &argv,
        MPI_THREAD_FUNNELED,
        &thread_support);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    fftw_mpi_init();

    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());

    int N = 512;
    string fname = "/store/uzh/uzh8/ESC412/ic_512.std";

    array2D_r p = read_particles(fname, rank, size);

    ptrdiff_t alloc_local, local_n0, local_0_start;

    // Compute local sizes
    alloc_local = fftw_mpi_local_size_3d(N, N, N, MPI_COMM_WORLD,
                                         &local_n0, &local_0_start);

    array2D_r p_balanced = balance_particles(p, local_0_start, N, rank, size);
    p.free();

    real_type* local_slab_real = fftw_alloc_real((local_n0+3)*N*2*(N/2+1));
    blitz::GeneralArrayStorage<3> storage;
    storage.base() = local_0_start, 0, 0;
    array3D_r grid(local_slab_real, blitz::shape(local_n0+3,N,2*(N/2+1)), storage);
    grid = 0;

    assign_masses(4, p_balanced, grid, rank, size);
    p_balanced.free();

    // Allocate the output buffer for the fft
    fftw_complex* local_slab_cplx = fftw_alloc_complex((local_n0)*N*(N/2+1));
    array3D_c fft_grid(reinterpret_cast<complex_type*>(local_slab_cplx),
                       blitz::shape(local_n0,N,N/2+1), storage);
    fft_grid = 0;

    // Compute the fft of the over-density field
    // The results are stored in fft_grid (out-of-place method)
    compute_fft(grid, fft_grid, N, MPI_COMM_WORLD);

    // Compute the power spectrum
    compute_pk(fft_grid, N, rank);

    fftw_mpi_cleanup();
    MPI_Finalize();

}
