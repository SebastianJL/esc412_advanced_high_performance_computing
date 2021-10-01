 // Compile: g++ -std=c++11 -O -Iinclude -o carray carray.cxx

#include <iostream>
#include "blitz/array.h"

int main() {
    constexpr int n = 30;

    typedef blitz::Array<int,1> array;

    array a(n);

    // set them to values 0 to n-1
    for(auto i=0; i<n; ++i) a(i) = i;

    for(auto i=0; i<n; ++i)
    	std::cout << i << " " << a(i) << std::endl;
    }