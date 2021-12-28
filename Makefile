mpi_openmp: src/main_mpi_openmp.cpp
	mpiCC -std=c++17 -I include/ -L /usr/lib/x86_64-linux-gnu/ -O2 src/main_mpi_openmp.cpp src/tipsy.cpp -o make-build/main_mpi_openmp -lfftw3_mpi -lfftw3_omp -lfftw3 -lm -lstdc++ -fopenmp

mpi_openmp_sample: src/main_mpi_openmp.cpp
	mpiCC -std=c++17 -I include/ -L /usr/lib/x86_64-linux-gnu/ -O2 src/sample_solution_mpi_openmp.cpp src/tipsy.cpp -o make-build/main_mpi_openmp_sample -lfftw3_mpi -lfftw3_omp -lfftw3 -lm -lstdc++ -fopenmp

cuda: src/main_cuda.cpp
	nvcc -std=c++17 -I include/  -O2 src/main_cuda.cpp src/tipsy.cpp -o make-build/main_cuda -lcudart -lcufft -lcufftw -lm -lstdc++

