#ifndef LLPS_CALCULUS_FINITE_DIFFERENCE_HPP_INCLUDED
#define LLPS_CALCULUS_FINITE_DIFFERENCE_HPP_INCLUDED

#include <array>    //Access to std::array
#include <utility>  //Access to std::swap
#include <iterator> //Access to std::iter_value_t and iterator concepts
#include <concepts> //Access to std::signed_integral concept
#include <vector>

namespace llps::calculus {

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
    constexpr auto fd_stencil(size_t order, std::array<IntType, N> samples)
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

            //Swapping i-th element with the last to emulate removing it. Note that range below is decreased by one.
            std::swap(samples[i], samples.back());
            result[i] = fact * elem_sym_poly(N - order - 1, samples.begin(), std::prev(samples.end())) / denom;
            //Revert samples back to unaltered state
            std::swap(samples[i], samples.back());
        }

        return result;
    }


    template<size_t error_order>
    consteval auto central_indicies()
    {
        static_assert((error_order & 1) == 0, "error order must be even for central finite difference!");

        std::array<ptrdiff_t, error_order + 1> result{};

        for (size_t i = 0; i <= error_order; ++i)
            result[i] = i - static_cast<ptrdiff_t>(error_order)/2;

        return result;
    }

    template<size_t error_order, class OutType = double>
    constexpr auto central_fd_stencil(size_t order)
    {
        return fd_stencil<OutType>(order, central_indicies<error_order>());
    }

}

#endif // !CALCULUS_FINITE_DIFFERENCE_HPP_INCLUDED