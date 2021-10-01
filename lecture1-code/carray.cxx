// Compile: g++ -std=c++11 -O -o carray carray.cxx

#include <iostream>

int main() {
    constexpr int n = 30;

    // allocate "n" integers
    // int * a = malloc(n * sizeof(int));
    // int * a = new int[n];
    auto a = new int[n];

    // set them to values 0 to n-1
    for(auto i=0; i<n; ++i) a[i] = i;

    for(auto i=0; i<n; ++i)
    	std::cout << i << " " << a[i] << std::endl;
    }