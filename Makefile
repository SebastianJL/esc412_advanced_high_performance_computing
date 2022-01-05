mpi_openmp: src/main_mpi_openmp.cpp src/tipsy.cpp
	mpiCC -std=c++17 -I include/ -L /usr/lib/x86_64-linux-gnu/ -O2 $^ -o make-build/$@ -lfftw3_mpi -lfftw3_omp -lfftw3 -lm -lstdc++ -fopenmp

mpi_openmp_sample: src/main_mpi_openmp_sample.cpp src/tipsy.cpp
	mpiCC -std=c++17 -I include/ -L /usr/lib/x86_64-linux-gnu/ -O2 $^ -o make-build/$@ -lfftw3_mpi -lfftw3_omp -lfftw3 -lm -lstdc++ -fopenmp

clean:
	rm -f make-build/*
