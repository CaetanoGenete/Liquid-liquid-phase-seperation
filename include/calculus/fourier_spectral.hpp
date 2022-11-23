#ifndef LLPS_CALCULUS_FOURIER_SPECTRAL_HPP_INCLUDED
#define LLPS_CALCULUS_FOURIER_SPECTRAL_HPP_INCLUDED

#include <array>
#include <cstddef>

namespace llps::calculus {

    template<size_t _rows>
    consteval std::array<uint32_t, _rows> row_freq_indicies()
    {
        std::array<uint32_t, _rows> result{};
        result[_rows / 2] = _rows / 2;

        for (size_t i = 0; i < _rows / 2; ++i) {
            result[i] = i;
            result[_rows - 1 - i] = i + 1;
        }

        return result;
    }

}

#endif // !LLPS_CALCULUS_FOURIER_SPECTRAL_HPP_INCLUDED