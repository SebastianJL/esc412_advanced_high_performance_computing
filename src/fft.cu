#include <fstream>
#include <typeinfo>
#include <complex>
#include <stdlib.h>
#include "cufft.h"
#include "blitz/array.h"

#include "aweights.h"
#include "tipsy.h"

using namespace std;
using namespace blitz;

typedef double real_type;
typedef std::complex<real_type> complex_type;
typedef blitz::Array<real_type,3> array3D_r;
typedef blitz::Array<complex_type,3> array3D_c;


//**********************************************************************

void compute_fft_2D_R2C(array3D_r &grid, int N, int local_n) {
    int n[] = {N,N};       // 2D FFT of length NxN
    int inembed[] = {N,2*(N/2+1)};
    int onembed[] = {N,(N/2+1)};
    int howmany = local_n;
    int odist = N*(N/2+1); // Output distance is in "complex"
    int idist = 2*odist;   // Input distance is in "real"
    int istride = 1;       // Elements of each FFT are adjacent
	int ostride = 1;

    cufftHandle plan;
    cufftPlanMany(&plan,sizeof(n)/sizeof(n[0]), n,
		    inembed,istride,idist,
		    onembed,ostride,odist,
		    CUFFT_D2Z,howmany);
    cufftDoubleComplex *data;
    auto data_size = sizeof(cufftDoubleComplex)*local_n*N*(N/2+1);
    cudaMalloc((void**)&data, data_size);
    cudaMemcpy(data, grid.dataFirst(), data_size, cudaMemcpyHostToDevice);
    cufftExecD2Z(plan,reinterpret_cast<cufftDoubleReal*>(data),data);
    cudaMemcpy(grid.dataFirst(), data, data_size, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    cudaFree(data);
    cufftDestroy(plan);

}

/** 1D complex to complex fft along the second dimension.
 *
 * In the complex-space the grid is N*N*(N/2+1). Since the code is parallelized
 * along the first dimension the actual grid size is local_n*N*(N/2+1).
 * It is assumed that compute_fft_2D_R2C and transpose() have already been
 * called on the grid.
 *
 * @param fft_grid
 * @param N Grid parameter.
 * @param local_n Local size of the first dimension.
 */
void compute_fft_1D_C2C(array3D_c &fft_grid, int N, int local_n){
    int n[] = {N};
    int *inembed = n;  // Effectively ignored. See doc for outermost dimension.
    int *onembed = n;  // Effectively ignored.
    int howmany = N;
    int idist = 1;
    int odist = 1;
    int istride = (N/2+1);
    int ostride = (N/2+1);

    
    for (int slab=0; slab<local_n; slab++) {

        cufftHandle plan;
        cufftPlanMany(&plan,sizeof(n)/sizeof(n[0]), n,
                      inembed,istride,idist,
                      onembed,ostride,odist,
                      CUFFT_Z2Z,howmany);
        
        if (cudaGetLastError() != cudaSuccess){
            fprintf(stderr, "Cuda error: Failed to plan\n");
        }

        cufftDoubleComplex *data;
        auto data_size = sizeof(cufftDoubleComplex)*N*(N/2+1);
        cudaMalloc((void**)&data, data_size);

        if (cudaGetLastError() != cudaSuccess){
            fprintf(stderr, "Cuda error: Failed to allocate\n");
        }

        int index = slab * N*(N/2 + 1);

        cudaMemcpy(data, fft_grid.dataFirst() + index, data_size, cudaMemcpyHostToDevice);
        if (cudaGetLastError() != cudaSuccess){
            fprintf(stderr, "Cuda error: Failed to copy to device\n");
        }
        
        cufftExecZ2Z(plan, data, data, CUFFT_FORWARD);
        
        if (cudaGetLastError() != cudaSuccess){
            fprintf(stderr, "Cuda error: Failed to execute plan\n");
        }

        cudaMemcpy(fft_grid.dataFirst() + index, data, data_size, cudaMemcpyDeviceToHost);
        if (cudaGetLastError() != cudaSuccess){
            fprintf(stderr, "Cuda error: Failed to copyt to host\n");
        }

        cudaDeviceSynchronize();
        cudaFree(data);
        cufftDestroy(plan);
    }

}

//void compute_fft_1D_C2C(array3D_c &fft_grid, int N, int local_n){
//    int n[] = {N};
//    int *inembed = n, *onembed = n;
//    int howmany = N*(N/2+1);
//    int idist = 1;
//    int odist = 1;
//    int istride = N*(N/2+1);
//    int ostride = N*(N/2+1);
//
//    cufftHandle plan;
//    cufftPlanMany(&plan,sizeof(n)/sizeof(n[0]), n,
//                    inembed,istride,idist,
//                    onembed,ostride,odist,
//                    CUFFT_Z2Z,howmany);
//    cufftDoubleComplex *data;
//    cudaMalloc((void**)&data, sizeof(cufftDoubleComplex)*N*N*(N/2+1));
//    cudaMemcpy(data, fft_grid.dataFirst(), sizeof(cufftDoubleComplex)*N*N*(N/2+1), cudaMemcpyHostToDevice);
//    cufftExecZ2Z(plan,data,data,CUFFT_FORWARD);
//    cudaMemcpy(fft_grid.dataFirst(), data,sizeof(cufftDoubleComplex)*N*N*(N/2+1), cudaMemcpyDeviceToHost);
//    cudaDeviceSynchronize();
//    cudaFree(data);
//    cufftDestroy(plan);
//
//}

//**********************************************************************

void compute_fft_2D_R2C_stream(array3D_r &grid, array3D_c &fft_grid, int N) {
    int n[] = {N,N};       // 2D FFT of length NxN
    int inembed[] = {N,2*(N/2+1)};
    int onembed[] = {N,(N/2+1)};
    const int howmany = 16;// Number of slabs to do at once
    int odist = N*(N/2+1); // Output distance is in "complex"
    int idist = 2*odist;   // Input distance is in "real"
    int istride = 1,       // Elements of each FFT are adjacent
	ostride = 1;
    const int nStreams = 4;

    // Allocate the CUDA streams. Each stream can execute independently
    cudaStream_t stream[nStreams];
    for(auto i=0; i<nStreams; ++i) cudaStreamCreate(stream+i);

    // Allocate a single chunk on the GPU, but separate it into blocks
    cufftDoubleComplex *data[nStreams];
    int block_count = howmany*N*(N/2+1);
    cudaMalloc((void**)data, sizeof(cufftDoubleComplex)*nStreams*block_count);
    for(auto i=1; i<nStreams; ++i) data[i] = data[0] + i*block_count;

    // Create a plan to do "howmany" slabs at a time
    // This plan will be run simultaneously on multiple streams so that
    // means we have to create "work areas" for each. Instead of calling
    // cufftPlanMany we call cufftMakePlanMany so we can disable auto allocation
    // of the work area. Later we need to call cufftSetWorkArea.
    cufftHandle plan;
    size_t workSize;
    cufftCreate(&plan);
    cufftSetAutoAllocation(plan,0);
    cufftMakePlanMany(plan,sizeof(n)/sizeof(n[0]), n,
		    inembed,istride,idist,
		    onembed,ostride,odist,
		    CUFFT_D2Z,howmany,&workSize);
    void *workArea[nStreams];
    // We allocate "nStreams" work areas and set pointer to them for each stream
    cudaMalloc(&workArea[0],workSize*nStreams);
    for(auto i=1; i<nStreams; ++i) workArea[i] = reinterpret_cast<char*>(workArea[0]) + i*workSize;


    // Distribute the work on the streams
    int iStream = 0;
    for(auto i=grid.lbound(0); i<=grid.ubound(0); i+=howmany) {
	cudaMemcpyAsync(data[iStream], &grid(i,0,0), sizeof(cufftDoubleComplex)*block_count, cudaMemcpyHostToDevice,stream[iStream]);
	cufftSetStream(plan,stream[iStream]);
	cufftSetWorkArea(plan,workArea[iStream]);
	cufftExecD2Z(plan,reinterpret_cast<cufftDoubleReal*>(data[iStream]),data[iStream]);
	cudaMemcpyAsync(&grid(i,0,0),data[iStream],sizeof(cufftDoubleComplex)*block_count, cudaMemcpyDeviceToHost,stream[iStream]);
	if (++iStream == nStreams) iStream = 0;
    }
    cudaDeviceSynchronize(); // Wait for all streams to complete
    cudaFree(data[0]);
    cudaFree(workArea[0]);
    cufftDestroy(plan);
    for(auto i=0; i<nStreams; ++i) cudaStreamDestroy(stream[i]);
}



//**********************************************************************
