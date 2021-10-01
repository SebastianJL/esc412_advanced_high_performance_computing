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

    for(auto i=a.begin(); i!=a.end(); ++i)
        std::cout << i.position() << " " << *i << std::endl;
    }
