#include "fftw/fftw3.h"

#include <numeric>
#include <iostream>
#include <numbers>
#include <algorithm>

#include "grid.hpp"

int main()
{
    using std::numbers::pi;

    static constexpr size_t n = 16;
    static constexpr double L = 2. * pi;

    double* nw_sqr = fftw_alloc_real(n/2 + 1);

    for (ptrdiff_t i = 0; i <= n/2; ++i)
    {
        double w = (2 * pi * i)/L;
        nw_sqr[i] = -(w*w);
    }

    double* in = nullptr;
    fftw_complex* out = nullptr;

    fftw_plan plan;

    in = fftw_alloc_real(n);
    out = fftw_alloc_complex(n/2 + 1);

    plan = fftw_plan_dft_r2c_1d(n, in, out, FFTW_PATIENT);
    
    for (size_t i = 0; i < n; ++i) {
        in[i] = std::exp(std::cos((i * L) / n));
        std::cout << in[i] << ", ";
    }
    std::cout << "\n_________________________\n";

    fftw_execute(plan);

    for (size_t i = 0; i < n / 2 + 1; ++i)
    {
        out[i][0] *= nw_sqr[i];
        out[i][1] *= nw_sqr[i];
    }

    fftw_plan back_plan;
    back_plan = fftw_plan_dft_c2r_1d(n, out, in, FFTW_PATIENT);

    fftw_execute(back_plan);

    for (size_t i = 0; i < n; ++i)
        std::cout << in[i]/n << ", ";
    std::cout << "\n_________________________\n";

    for (size_t i = 0; i < n; ++i) {
        double x = (i * L) / n;

        std::cout << std::exp(std::cos(x)) * (std::sin(x) * std::sin(x) - std::cos(x)) << ", ";
    }

    fftw_destroy_plan(plan);
    fftw_destroy_plan(back_plan);

    fftw_free(nw_sqr);
    fftw_free(in);
    fftw_free(out);
}