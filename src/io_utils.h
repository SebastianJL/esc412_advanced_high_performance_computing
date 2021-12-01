//
// Created by johannes on 11/4/21.
//

#ifndef PROJECT_IO_UTILS_H
#define PROJECT_IO_UTILS_H

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

template<typename t>
void print_vector(std::vector<t>& vec){
    std::cout << '[';
    for (int i = 0; i < vec.size(); i++)
        std::cout << vec[i] << " ";
    std::cout << ']' << std::endl;
}

template<typename t>
std::string sprint_array(t* arr, int length){
    std::stringstream buffer;
    buffer << "{ ";
    for (int i = 0; i < length; i++)
        buffer << arr[i] << " ";
    buffer << '}' << std::flush;
    return buffer.str();
}

#endif // PROJECT_IO_UTILS_H
