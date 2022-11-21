#ifndef CALCULUS_FINITE_DIFFERENCE_HPP_INCLUDED
#define CALCULUS_FINITE_DIFFERENCE_HPP_INCLUDED

#include <array>
#include <cstddef>
#include <numeric>
#include <iterator>

namespace calculus {

    template<std::forward_iterator It>
    constexpr auto elem_sym_poly(size_t k, It first, It last)
    {
        const size_t N = std::distance(first, last);
        const size_t M = N - k + 1;

        //Todo: perhaps avoid need for 1 here
        std::vector<std::iter_value_t<It>> e_n(M, 1);

        for (size_t i = 0; i < k; ++i)
        {
            e_n[0] = *first * e_n[0];

            auto it = ++first;
            for (size_t j = 1; j < M; ++j, ++it)
                e_n[j] = *it * e_n[j] + e_n[j - 1];
        }

        return e_n.back();
    }

    template<class OutType = double, std::signed_integral IntType, size_t N>
    constexpr auto finite_difference_stencil(size_t order, std::array<IntType, N> samples)
    {
        std::array<OutType, N> result;

        OutType fact((N - order) % 2 == 0 ? -1 : 1);
        for (size_t i = 2; i < order + 1; ++i)
            fact *= i;

        for (size_t i = 0; i < N; ++i)
        {
            OutType denom = 1;
            for (size_t j = 0; j < N; ++j) {
                if (j != i)
                    denom *= (samples[i] - samples[j]);
            }

            //Swapping here to emulate removing an element
            std::swap(samples[i], samples.back());
            result[i] = fact * elem_sym_poly(N - order - 1, samples.begin(), std::prev(samples.end())) / denom;
            //Revert samples back to unaltered state
            std::swap(samples[i], samples.back());
        }

        return result;
    }
}

#endif // !CALCULUS_FINITE_DIFFERENCE_HPP_INCLUDED