//
// Created by johannes on 1/19/22.
//

#include "get_time.h"

// A simple timing function
double get_time() {
    struct timeval tv;
    struct timezone tz;
    gettimeofday(&tv, &tz);
    return tv.tv_sec + 1e-6*(double)tv.tv_usec;
}
