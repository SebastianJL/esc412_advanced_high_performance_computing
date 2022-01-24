# esc412_advanced_high_performance_computing


This project takes positional data of particles distributed similarly to star-clusters in the universe and computes the
power spectrum.

## Project structure
Code written by me can be found in src/.

The project contains multiple versions of the same code. There is a serial version (not finished) and several versions
that where parallelized using OpenMP, MPI and Cuda. This is refelected in the name of the main files:
- main_serial.cpp
- main_mpi_openmp.cpp
- main_mpi_openmp_sample.cpp (A file containing sample code given by the assistants.)
- main_cuda.cu
- main_cuda_mpi_openmp.cpp


## Build the project

- First build blitz.
    ```bash
    cmake -DCMAKE_INSTALL_PREFIX=./ -B ./blitz-build ./blitz-1.0.2
    make -C ./blitz-build
    make -C ./blitz-build install
    ```
    This will create the folders `include`, `lib` and `blitz-build`. They need to be present in order for the
    project to run.

- The project itself can be built with one of the targets in `Makefile` or `Makefile.daint`. The binaries are built into
  `make-build`.
  Slurm scripts to run the compiled code on daint (or eiger) can be found in scripts/. They are named
  `run_<name_of_executable>`. 
- The cmake stuff can largely be ignored as it does not work and is only present such that the clion code introspection
  features work.