#include <cstdint> 
#include "tipsy.h"

int main() {
    using std::cout, std::endl;
    using namespace blitz;

    TipsyIO io;

    io.open("b0-final.std");
    std::cout << "N = " << io.count() << endl << endl;

    if (io.fail()) {
        std::cerr << "Unable to open file" << std::endl;
        abort();
    }

    Array<int, 2> A(4,4);
    A = 3, 8, 0, 1,
        1, -1, 9, 3,
        2, -5, -1, 1,
        4, 3, 4, 2;


    Array<float,2> r(io.count(),3);
    io.load(r);
    
    cout << r(Range(0,10), Range::all()) << endl;

    secondIndex j;
    blitz::Array<float, 1> m = blitz::max(r, j);
    cout << "max:\n" << m << endl;
}
