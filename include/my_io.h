//
// Created by johannes on 11/4/21.
//

#ifndef PROJECT_MY_IO_H
#define PROJECT_MY_IO_H

#include "blitz/array.h"
template <typename real_t>
int write_to_csv(blitz::Array<real_t, 2> array, const char *filename);

#endif // PROJECT_MY_IO_H
