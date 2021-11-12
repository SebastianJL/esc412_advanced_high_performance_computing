//
// Created by johannes on 11/4/21.
//

#ifndef PROJECT_MY_IO_H
#define PROJECT_MY_IO_H

#include "blitz/array.h"
template <typename real_t>
int write_to_csv(blitz::Array<real_t, 2> array, const char *filename);

template<typename real_t>
blitz::Array<real_t, 2>
project_3d_to_2d(const blitz::Array<real_t, 3> &grid_3d) {
    int n_grid = grid_3d.extent(blitz::firstDim);
    blitz::thirdIndex k;
    blitz::Array<real_t, 2> grid_2d(n_grid, n_grid);

    grid_2d = max(grid_3d, k);
    return grid_2d;
}


#endif // PROJECT_MY_IO_H
