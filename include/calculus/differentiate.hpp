#ifndef CALCULUS_DIFFERENTIATE_HPP_INCLUDED
#define CALCULUS_DIFFERENTIATE_HPP_INCLUDED

#include <concepts> // For access to iterator concepts
#include <array>
#include <type_traits>

#include "finite_difference.hpp"
#include "fourier_spectral.hpp"
#include "../grid.hpp"

namespace calculus {
    /*
    template<typename Type, size_t _rows, size_t _cols, typename Container, std::integral IntType, size_t _kernel_size>
    constexpr auto laplacian_fd(
        const grid<Type, _rows, _cols, Container>& phi, 
        const std::array<IntType, _kernel_size>& samples,
        Type dx, Type dy)
    {
        static auto stencil = fd_stencil<Type>(2, samples);
        static constexpr size_t offset = _kernel_size / 2;

        grid<Type, _rows, _cols, Container> dphi;

        for (size_t row = 0; row < _rows; ++row)
        {
            size_t offset_row = (row + offset) % _rows;;

            for (size_t col = 0; col < _cols; ++col)
            {
                size_t offset_col = (col + offset) % _cols;

                for (size_t i = 0; i < _kernel_size; ++i) {
                    //Using 2D index here to remain as general as possible
                    dphi(offset_row, offset_col) += (
                        phi(offset_row, (col + i) % _cols) / (dx * dx) +
                        phi((row + i) % _rows, offset_col) / (dy * dy)
                    ) * stencil[i];
                }
            }
        }

        return dphi;
    }
    */

    template<size_t error_order, typename Type, size_t _rows, size_t _cols, typename Container>
    constexpr auto laplacian_central_fd(const grid<Type, _rows, _cols, Container>& phi, Type dx, Type dy)
    {
        static constexpr auto stencil = central_fd_stencil<error_order, Type>(2);
        static constexpr size_t offset = error_order / 2;

        grid<Type, _rows, _cols, Container> dphi;

        for (size_t row = 0; row < _rows; ++row)
        {
            size_t offset_row = (row + offset) % _rows;;

            for (size_t col = 0; col < _cols; ++col)
            {
                size_t offset_col = (col + offset) % _cols;

                for (size_t i = 0; i <= error_order; ++i) {
                    //Using 2D index here to remain as general as possible
                    dphi(offset_row, offset_col) += (
                        phi(offset_row, (col + i) % _cols) / (dx * dx) +
                        phi((row + i) % _rows, offset_col) / (dy * dy)
                    ) * stencil[i];
                }
            }
        }

        return dphi;
    }
}

#ifdef LLPS_USE_MKL

#include "fftw/fftw3.h"

namespace calculus {

    template<typename Type>
    struct as_ftw_complex;

    template<> 
    struct as_ftw_complex<float> { using value_type = fftwf_complex; };
    template<>
    struct as_ftw_complex<double> { using value_type = fftw_complex; };
    template<>
    struct as_ftw_complex<long double> { using value_type = fftwl_complex; };

    template<typename Type>
    using as_ftw_complex_t = typename as_ftw_complex<Type>::value_type;

    template<std::floating_point Type, size_t _rows, size_t _cols, class Container1, class Container2>
    void laplacian_spectral(
        grid<Type, _rows, _cols, Container1>& phi,
        grid<Type, _rows, _cols, Container2>& dphi,
        Type dx, Type dy)
    {
        using complex_type = as_ftw_complex_t<Type>;

        static constexpr size_t phi_hat_size = _rows * (_cols/2 + 1);
        complex_type* phi_hat = static_cast<complex_type*>(fftw_malloc(sizeof(complex_type) * phi_hat_size));

        Type x_freq_elem = 2. * std::numbers::pi / (_cols * dx);
        Type y_freq_elem = 2. * std::numbers::pi / (_rows * dy);

        fftw_plan forw_plan = fftw_plan_dft_r2c_2d(_rows, _cols, phi.data(), phi_hat, FFTW_PATIENT);
        fftw_plan back_plan = fftw_plan_dft_c2r_2d(_rows, _cols, phi_hat, dphi.data(), FFTW_PATIENT);
        fftw_execute(forw_plan);

        static constexpr auto row_indicies = calculus::row_freq_indicies<_rows>();

        size_t index = 0;
        for (size_t row : row_indicies)
        {
            Type kappa_y = y_freq_elem * row;

            for (size_t col = 0; col <= _cols / 2; ++col, ++index)
            {
                Type kappa_x = x_freq_elem * col;
                Type kappa_xy = (-kappa_x * kappa_x - kappa_y * kappa_y)/ (_rows * _cols);

                phi_hat[index][0] *= kappa_xy;
                phi_hat[index][1] *= kappa_xy;
            }
        }

        fftw_execute(back_plan);

        fftw_destroy_plan(forw_plan);
        fftw_destroy_plan(back_plan);
        fftw_free(phi_hat);
    }
    template<std::floating_point Type, size_t _rows, size_t _cols, class Container>
    void laplacian_spectral(grid<Type, _rows, _cols, Container>& phi, Type dx, Type dy)
    {
        laplacian_spectral(phi, phi, dx, dy);
    }

}

#endif // LLPS_USE_MKL

#endif // !CALCULUS_DIFFERENTIATE_HPP_INCLUDED