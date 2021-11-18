//
// Created by johannes on 11/13/21.
//

#ifndef PROJECT_FOURIER_TRANSFORM_H
#define PROJECT_FOURIER_TRANSFORM_H

#include "fftw3.h"

template <typename real_t, typename complex_t>
void compute_fft(blitz::Array<real_t, 3> grid,
                 blitz::Array<complex_t, 3> fft_grid, int n_grid) {
    auto plan = fftw_plan_dft_r2c_3d(
        n_grid, n_grid, n_grid, grid.dataFirst(),
        reinterpret_cast<fftw_complex *>(fft_grid.dataFirst()), FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
}
#endif // PROJECT_FOURIER_TRANSFORM_H
