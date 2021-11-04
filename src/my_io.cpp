//
// Created by johannes on 11/4/21.
//

#include "my_io.h"
#include "blitz/array.h"
#include <fstream>

template <typename real_t>
int write_to_csv(blitz::Array<real_t, 2> array, const char *filename) {
    std::ofstream out;
    out.open(filename);
    if (!out.is_open()) {
        return 1;
    }

    int i_max = array.extent(blitz::firstDim);
    int j_max = array.extent(blitz::secondDim);
    for (auto i = 0; i < i_max; ++i) {
        for (auto j = 0; j < j_max - 1; ++j) {
            out << array(i, j) << ",";
        }
        out << array(i, j_max - 1);
        out << std::endl;
    }
    out.close();
    return 0;
}

// Declare concrete types of template.
template int write_to_csv(blitz::Array<double, 2> array, const char *filename);
template int write_to_csv(blitz::Array<float, 2> array, const char *filename);
