#include "gtest/gtest.h"

#include <iostream>
#include <tuple>
#include <algorithm>

#include "calculus/finite_difference.hpp"
#include "utilities/meta.hpp"

TEST(finite_difference_tests, test_equidistant_central_difference_stencil)
{
    //Data taken from: https://en.wikipedia.org/wiki/Finite_difference_coefficient

    static constexpr auto expected_matrix = std::make_tuple(
        std::make_tuple(
            std::array{-1./2., 0., 1./2.},
            std::array{1., -2., 1.}
        ),
        std::make_tuple(
            std::array{1./12., -2./3., 0., 2./3., -1./12.},
            std::array{-1./12., 4./3., -5./2., 4./3., -1./12.},
            std::array{-1./2., 1., 0., -1., 1./2.},
            std::array{1., -4., 6., -4., 1.}
        ),
        std::make_tuple(
            std::array{-1./60., 3./20., -3./4., 0., 3./4., -3./20., 1./60.},
            std::array{1./90., -3./20., 3./2., -49./18., 3./2., -3./20., 1./90.},
            std::array{1./8., -1., 13./8., 0., -13./8., 1., -1./8.},
            std::array{-1./6., 2., -13./2., 28./3., -13./2., 2., -1./6.},
            std::array{-1./2., 2., -5./2., 0., 5./2., -2., 1./2.},
            std::array{1., -6., 15., -20., 15., -6., 1.}
        ),
        std::make_tuple(
            std::array{1./280., -4./105., 1./5., -4./5., 0., 4./5., -1./5., 4./105., -1./280.},
            std::array{-1./560., 8./315., -1./5., 8./5., -205./72., 8./5., -1./5., 8./315., -1./560.},
            std::array{-7./240., 3./10., -169./120., 61./30., 0., -61./30., 169./120., -3./10., 7./240.},
            std::array{7./240., -2./5., 169./60., -122./15., 91./8., -122./15., 169./60., -2./5., 7./240.},
            std::array{1./6.,  -3./2., 13./3., -29./6., 0., 29./6., -13./3., 3./2., -1./6.},
            std::array{-1./4., 3., -13., 29., -75./2., 29., -13., 3., -1./4.}
        )
    );

    static constexpr size_t indicies_count = std::tuple_size_v<decltype(expected_matrix)>;

    llps::utilities::constexpr_for<indicies_count>([]<size_t I>(llps::utilities::size_t_constant<I>) {
        constexpr size_t max_order = std::min<size_t>(I * 2 + 2, 6);

        llps::utilities::constexpr_for<max_order>([]<size_t order>(llps::utilities::size_t_constant<order>) {
            auto& expected = std::get<order>(std::get<I>(expected_matrix));
            //Compiler just gives up if I set this to constexpr, LOL.
            auto actual = llps::calculus::central_fd_stencil<(I+1)*2>(order + 1);


            ASSERT_TRUE(std::ranges::equal(expected, actual)) << "Failed at: I=" << I << ", order=" << order;
        });
    });
}

TEST(finite_difference_tests, test_equidistant_forward_difference_stencil)
{
    //Data taken from: https://en.wikipedia.org/wiki/Finite_difference_coefficient

    static constexpr auto indicies = std::make_tuple(
        std::array{0, 1},
        std::array{0, 1, 2},
        std::array{0, 1, 2, 3},
        std::array{0, 1, 2, 3, 4});

    static constexpr size_t indicies_size = std::tuple_size_v<decltype(indicies)>;

    static constexpr auto expected_matrix = std::make_tuple(
        std::make_tuple(
            std::array{ -1., 1. }
        ),
        std::make_tuple(
            std::array{ -3./2., 2., -1./2. },
            std::array{ 1., -2., 1.}
        ),
        std::make_tuple(
            std::array{ -11./6., 3., -3./2., 1./3. },
            std::array{ 2., -5., 4., -1.},
            std::array{ -1., 3., -3., 1.}
        ),
        std::make_tuple(
            std::array{ -25./12., 4., -3., 4./3., -1./4. },
            std::array{ 35./12., -26./3., 19./2., -14./3., 11./12.},
            std::array{ -5./2., 9., -12., 7., -3./2. },
            std::array{1., -4., 6., -4., 1.}
        )
    );

    llps::utilities::constexpr_for<indicies_size>([]<size_t I>(llps::utilities::size_t_constant<I>) {
        llps::utilities::constexpr_for<I+1>([]<size_t order>(llps::utilities::size_t_constant<order>) {
            auto& expected = std::get<order>(std::get<I>(expected_matrix));
            //Compiler just gives up if I set this to constexpr, LOL.
            auto actual = llps::calculus::fd_stencil(order + 1, std::get<I>(indicies));

            ASSERT_TRUE(std::ranges::equal(expected, actual)) << "Failed at: I=" << I << ", order=" << order;
        });
    });
}


TEST(finite_difference_tests, test_central_indicies)
{
    static constexpr size_t test_count = 10;

    llps::utilities::constexpr_for<test_count>([]<size_t I>(llps::utilities::size_t_constant<I>) {
        constexpr size_t order = (I + 1) * 2;
        auto indicies = llps::calculus::central_indicies<order>();

        ASSERT_EQ(order + 1, std::size(indicies)) << "Mismatch in expected size!";

        for (size_t i = 0; i <= order; ++i)
            ASSERT_EQ(indicies[i], static_cast<ptrdiff_t>(i - order/2));
    });
}