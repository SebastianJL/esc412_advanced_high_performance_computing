// Compile: g++ -std=c++11 -O -Iinclude -o carray carray.cxx

#include <iostream>
#include "blitz/array.h"
using namespace blitz;
using namespace std;

int main() {
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

    //array::iterator i;
    array::const_iterator i;
    for(i=a.begin(); i!=a.end(); ++i)
        cout << i.position() << " " << *i << endl;
    }