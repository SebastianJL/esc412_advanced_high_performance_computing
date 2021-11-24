cmake_minimum_required(VERSION 3.19)
project(power_spectrum)

set(CMAKE_CXX_STANDARD 17)

include_directories(src)

add_executable(test_mass_assignement
        src/test_mass_assignement.cpp
        src/tipsy.cpp)
target_link_libraries(test_mass_assignement blitz)

add_executable(power_spectrum
        src/main.cpp
        src/tipsy.cpp
        src/mass_assignement.h
        src/fourier_transform.h
        src/io_utils.h)
target_link_libraries(power_spectrum blitz)
target_link_libraries(power_spectrum fftw3)


add_executable(main_ex9
        src/main_ex9.cpp
        src/tipsy.cpp
        )
target_link_libraries(main_ex9 blitz)
target_link_libraries(main_ex9 fftw3-mpi)