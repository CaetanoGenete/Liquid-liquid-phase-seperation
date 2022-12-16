#ifndef LLPS_CALCULUS_DIFFERENTIATE_HPP_INCLUDED
#define LLPS_CALCULUS_DIFFERENTIATE_HPP_INCLUDED

#include <concepts> // For access to iterator concepts
#include <array>
#include <type_traits>
#include <numbers>

#include "finite_difference.hpp"
#include "fourier_spectral.hpp"
#include "../grid.hpp"

namespace llps::calculus {

    template<size_t error_order, grid_like InGrid, grid_like OutGrid>
    LLPS_FORCE_INLINE constexpr void laplacian_central_fd(
        const InGrid& phi, OutGrid& dphi, 
        typename OutGrid::value_type dx, 
        typename OutGrid::value_type dy)
    {
        static constexpr auto stencil = central_fd_stencil<error_order, typename OutGrid::value_type>(2);
        static constexpr size_t offset = error_order / 2;

        for (size_t row = 0; row < phi.rows(); ++row)
        {
            size_t offset_row = (row + offset) % phi.rows();

            for (size_t col = 0; col < phi.cols(); ++col)
            {
                size_t offset_col = (col + offset) % phi.cols();

                dphi(offset_row, offset_col) = 0.;
                for (size_t i = 0; i <= error_order; ++i) {
                    //Using 2D index here to remain as general as possible
                    dphi(offset_row, offset_col) += (
                        phi(offset_row, (col + i) % phi.cols()) / (dx * dx) +
                        phi((row + i) % phi.rows(), offset_col) / (dy * dy)
                    ) * stencil[i];
                }
            }
        }
    }

    template<size_t error_order, class Meta>
    LLPS_FORCE_INLINE constexpr auto laplacian_central_fd(
        const llps::_basic_grid<Meta>& phi, 
        const llps::grid_value_t<Meta> dx,
        const llps::grid_value_t<Meta> dy)
    {
        llps::_basic_grid<Meta> dphi;
        laplacian_central_fd<error_order>(phi, dphi, dx, dy);

        return dphi;
    }
}

#ifdef LLPS_USE_MKL

#include "fftw/fftw3.h"

namespace llps::calculus {

    template<typename Type>
    struct as_ftw_complex;

    template<> struct as_ftw_complex<float> { using value_type = fftwf_complex; };
    template<> struct as_ftw_complex<double> { using value_type = fftw_complex; };
    template<> struct as_ftw_complex<long double> { using value_type = fftwl_complex; };

    template<typename Type>
    struct complex_base;

    template<> struct complex_base<fftwf_complex> { using value_type = float; };
    template<> struct complex_base<fftw_complex> { using value_type = double; };
    template<> struct complex_base<fftwl_complex> { using value_type = long double; };

    template<typename Type>
    using as_ftw_complex_t = typename as_ftw_complex<Type>::value_type;

    template<typename Type>
    using complex_base_t = typename complex_base<Type>::value_type;

    template<size_t _rows, size_t _cols, std::input_iterator It>
    void mult_herm_nfreq_squared(
        It first, 
        complex_base_t<std::iter_value_t<It>> dx,
        complex_base_t<std::iter_value_t<It>> dy)
    {
        using type = complex_base_t<std::iter_value_t<It>>;

        type x_freq_elem = 2. * std::numbers::pi / (_cols * dx);
        type y_freq_elem = 2. * std::numbers::pi / (_rows * dy);

        static constexpr auto row_indicies = row_freq_indicies<_rows>();

        for (int32_t row : row_indicies)
        {
            type kappa_y = y_freq_elem * row;

            for (size_t col = 0; col <= _cols / 2; ++col, ++first)
            {
                type kappa_x = x_freq_elem * col;
                type kappa_xy = (-kappa_x * kappa_x - kappa_y * kappa_y) / (_rows * _cols);

                (*first)[0] *= kappa_xy;
                (*first)[1] *= kappa_xy;
            }
        }
    }

    template<std::floating_point Type, size_t _rows, size_t _cols, class Container1, class Container2>
    void laplacian_spectral(
        llps::grid<Type, _rows, _cols, Container1>& phi,
        llps::grid<Type, _rows, _cols, Container2>& dphi,
        Type dx, Type dy)
    {
        using complex_type = as_ftw_complex_t<Type>;

        static constexpr size_t phi_hat_size = _rows * (_cols/2 + 1);
        complex_type* phi_hat = static_cast<complex_type*>(fftw_malloc(sizeof(complex_type) * phi_hat_size));

        fftw_plan forw_plan = fftw_plan_dft_r2c_2d(_rows, _cols, phi.data(), phi_hat, FFTW_ESTIMATE);
        fftw_plan back_plan = fftw_plan_dft_c2r_2d(_rows, _cols, phi_hat, dphi.data(), FFTW_ESTIMATE);
        
        fftw_execute(forw_plan);
        mult_herm_nfreq_squared<_rows, _cols>(phi_hat, dx, dy);
        fftw_execute(back_plan);

        fftw_destroy_plan(forw_plan);
        fftw_destroy_plan(back_plan);
        fftw_free(phi_hat);
    }

    template<std::floating_point Type, size_t _rows, size_t _cols, class Container>
    void laplacian_spectral(llps::grid<Type, _rows, _cols, Container>& phi, Type dx, Type dy)
    {
        laplacian_spectral(phi, phi, dx, dy);
    }

}

#endif // LLPS_USE_MKL

#endif // !LLPS_CALCULUS_DIFFERENTIATE_HPP_INCLUDED