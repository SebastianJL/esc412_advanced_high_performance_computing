// Compile: g++ -std=c++11 -O -Iinclude -o carray carray.cxx

#include <iostream>
#include "blitz/array.h"

int main() {
    constexpr int nx = 6; // "rows"
    constexpr int ny = 5; // "columns"
    constexpr int n = nx * ny;

    // allocate "n" integers
    // int * a = new int[n];
    auto data = new int[n];

    // set them to values 0 to n-1
    for(auto i=0; i<n; ++i) data[i] = i;

    typedef blitz::Array<int,2> array;

    array a(data,blitz::shape(nx,ny),blitz::neverDeleteData);

    for(auto i=0; i<nx; ++i)
	for(auto j=0; j<ny; ++j)
            std::cout << i << " " << j << " " << a(i,j) << std::endl;
    }