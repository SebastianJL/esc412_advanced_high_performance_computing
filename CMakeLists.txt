cmake_minimum_required(VERSION 3.19)
project(power_spectrum)
include(CMakePrintHelpers)

set(CMAKE_CXX_STANDARD 17)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")


include_directories(include)
include_directories(include/blitz)
include_directories(include/blitz/array)
include_directories(include/blitz/meta)
include_directories(include/random)
include_directories(src)


# test_mass_assignement
add_executable(test_mass_assignement
        src/test_mass_assignement.cpp
        src/tipsy.cpp)


# power spectrum
add_executable(power_spectrum
        src/main.cpp
        src/tipsy.cpp
)
target_link_libraries(power_spectrum fftw3)


# MPI + OpenMP
find_package(MPI REQUIRED)
find_package(FFTW REQUIRED)
#set(MDL_FFTW ${FFTW_FOUND})
find_package(OpenMP REQUIRED)
add_executable(main_mpi_openmp
        src/main_mpi_openmp.cpp
        src/tipsy.cpp
        src/transformations.h
        )
cmake_print_variables(FFTW_INCLUDES)
cmake_print_variables(FFTW_LIBRARIES)
target_link_libraries(main_mpi_openmp PUBLIC ${MPI_C_LIBRARIES})
target_link_libraries(main_mpi_openmp PUBLIC ${FFTW_LIBRARIES})
if(OpenMP_CXX_FOUND)
    target_link_libraries(main_mpi_openmp PUBLIC OpenMP::OpenMP_CXX)
endif()
target_include_directories(main_mpi_openmp PUBLIC ${MPI_C_INCLUDE_PATH})
target_include_directories(main_mpi_openmp PUBLIC ${FFTW_INCLUDES})


# Test transformations.h
add_executable(test_transformations
        src/test_transformations.cpp
        )