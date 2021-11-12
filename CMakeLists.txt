cmake_minimum_required(VERSION 3.19)
project(project)

set(CMAKE_CXX_STANDARD 17)

include_directories(include)
include_directories(include/blitz)
include_directories(include/blitz/array)
include_directories(include/blitz/meta)
include_directories(include/random)
include_directories(src)

add_executable(test_mass_assignement
        src/test_mass_assignement.cpp
        src/tipsy.cpp)

add_executable(project
        src/main.cpp
        src/tipsy.cpp)
