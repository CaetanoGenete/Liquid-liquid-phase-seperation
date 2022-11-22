#ifndef CALCULUS_DIFFERENTIATE_HPP_INCLUDED
#define CALCULUS_DIFFERENTIATE_HPP_INCLUDED

#include <concepts> // For access to iterator concepts
#include <array>

#include "finite_difference.hpp"
#include "../grid.hpp"

namespace calculus {

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

    template<size_t error_order, typename Type, size_t _rows, size_t _cols, typename Container>
    constexpr auto laplacian_central_fd(const grid<Type, _rows, _cols, Container>& phi, Type dx, Type dy)
    {
        return laplacian_fd(phi, central_indicies<error_order>(), dx, dy);
    }

}

#endif // !CALCULUS_DIFFERENTIATE_HPP_INCLUDED