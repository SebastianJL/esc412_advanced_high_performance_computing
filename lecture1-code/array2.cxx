// Compile: g++ -std=c++11 -O -Iinclude -o carray carray.cxx

#include <iostream>
#include "blitz/array.h"

int main() {
    constexpr int nx = 6; // "rows"
    constexpr int ny = 5; // "columns"
    constexpr int n = nx * ny;

    typedef blitz::Array<int,2> array;

    array a(nx,ny);

    for(auto i=0; i<nx; ++i)
	for(auto j=0; j<ny; ++j)
            a(i,j) = i*10 + j;

    for(auto i=0; i<nx; ++i)
	for(auto j=0; j<ny; ++j)
            std::cout << i << " " << j << " " << a(i,j) << std::endl;
    }
