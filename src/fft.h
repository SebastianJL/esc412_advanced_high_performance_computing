//
// Created by johannes on 1/19/22.
//

#ifndef POWER_SPECTRUM_FFT_H
#define POWER_SPECTRUM_FFT_H

#include "blitz/array.h"

typedef double real_type;
typedef std::complex<real_type> complex_type;
typedef blitz::Array<real_type,3> array3D_r;
typedef blitz::Array<complex_type,3> array3D_c;

void compute_fft_2D_R2C(array3D_r &grid, int N, int local_n);
void compute_fft_1D_C2C(array3D_c &fft_grid, int N, int local_n);
void compute_fft_2D_R2C_stream(array3D_r &grid, array3D_c &fft_grid, int N);


#endif // POWER_SPECTRUM_FFT_H
