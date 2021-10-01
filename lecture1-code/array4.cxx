// Compile: g++ -std=c++17 -O -Iinclude -o carray carray.cxx

#include <iostream>
#include "blitz/array.h"

int main() {
    using std::cout, std::endl;
    using blitz::Array, blitz::shape, blitz::neverDeleteData, blitz::FortranArray;

    constexpr int nx = 6; // "rows"
    constexpr int ny = 5; // "columns"
    constexpr int n = nx * ny;

    // allocate "n" integers
    // int * a = new int[n];
    auto data = new int[n];

    // set them to values 0 to n-1
    for(auto i=0; i<n; ++i) data[i] = i;

    typedef Array<int,2> array;

    array a(data,shape(nx,ny),neverDeleteData,FortranArray<2>());

    for(auto i=1; i<=nx; ++i)
	for(auto j=1; j<=ny; ++j)
            cout << i << " " << j << " " << a(i,j) << endl;
    }