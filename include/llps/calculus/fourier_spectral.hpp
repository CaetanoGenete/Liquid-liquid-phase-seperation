#ifndef LLPS_CALCULUS_FOURIER_SPECTRAL_HPP_INCLUDED
#define LLPS_CALCULUS_FOURIER_SPECTRAL_HPP_INCLUDED

#include <array>
#include <cstddef>

namespace llps::calculus {
    template<size_t _rows>
    consteval std::array<int32_t, _rows> row_freq_indicies()
    {
        static_assert(_rows % 2 == 0, "_rows must be even.");

        std::array<int32_t, _rows> result;
        for (int32_t i = 0; i < _rows / 2; ++i) {
            result[i] = i;
            result[_rows - 1 - i] = -i-1;
        }

        return result;
    }

}

#endif // !LLPS_CALCULUS_FOURIER_SPECTRAL_HPP_INCLUDED