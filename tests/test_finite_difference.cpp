#include "gtest/gtest.h"

#include <iostream>
#include <tuple>
#include <algorithm>
#include <vector>    //Access to std::vector
#include <numbers>   //Access to std::numbers::pi
#include <ranges>    //Access to std::ranges::views::take

#include "calculus/finite_difference.hpp"
#include "utilities/data_analytics.hpp"
#include "calculus/differentiate.hpp"
#include "utilities/meta.hpp"
#include "grid.hpp"

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

//Todo: Move to EXPU

template<typename SeqIt, typename Proj = std::identity>
struct transform_sequence;

template<typename SeqIt, typename Proj = std::identity>
using transform_sequence_t = typename transform_sequence<SeqIt, Proj>::type;

template<typename IntType, IntType ... _ints, typename Proj>
struct transform_sequence<std::integer_sequence<IntType, _ints...>, Proj>
{
    using type = std::integer_sequence<std::invoke_result_t<Proj, IntType>, Proj{}(_ints)...>;
};


template<typename SeqIt>
struct int_seq_to_gtest_Types;

template<typename IntType, IntType ... _ints>
struct int_seq_to_gtest_Types<std::integer_sequence<IntType, _ints...>>
{
    using type = testing::Types<std::integral_constant<IntType, _ints>...>;
};

template<typename SeqIt>
using int_seq_to_gtest_Types_t = typename int_seq_to_gtest_Types<SeqIt>::type;

template<size_t N, size_t M = 1, size_t C = 0>
using make_linear_index_sequence = transform_sequence_t < 
    std::make_index_sequence<N>, decltype([](size_t i) {
        return M * i + C;
    })>;


struct _fd_convergence_tests_data 
{
    static constexpr size_t N = 7;
    static constexpr size_t M = 2;
    static constexpr size_t C = 2;
};

template<class IntConst>
struct fd_convergence_tests;

template<size_t _n>
struct fd_convergence_tests<std::integral_constant<size_t, _n>> : 
    public testing::Test, _fd_convergence_tests_data
{
    static constexpr size_t error_order = _n;
};

using fd_convergence_test_types = int_seq_to_gtest_Types_t<
    make_linear_index_sequence<
        _fd_convergence_tests_data::N, 
        _fd_convergence_tests_data::M,
        _fd_convergence_tests_data::C>>;

TYPED_TEST_SUITE(fd_convergence_tests, fd_convergence_test_types);


TYPED_TEST(fd_convergence_tests, test_convergence)
{
    constexpr size_t test_index = (TestFixture::error_order - TestFixture::C) / TestFixture::M;
    constexpr double expected = (std::array{ 1.98884, 3.96378, 5.91638, 7.68905, 9.0925, 10.842, 12.1835})[test_index] ;

    using value_type = double;
    struct test_func
    {
        static value_type phi(value_type x, value_type y)
        {
            return std::exp(std::cos(x) + std::sin(y));
        }

        static value_type dphi(value_type x, value_type y)
        {

            return phi(x, y) * (std::cos(y) * std::cos(y) + std::sin(x) * std::sin(x) - std::sin(y) - std::cos(x));
        }
    };

    static constexpr size_t samples = 30;
    static constexpr value_type x_min = 0;
    static constexpr value_type x_max = 2. * std::numbers::pi;

    std::vector<value_type> delta_xs;
    std::vector<value_type> max_abs_errs;
    delta_xs.reserve(samples);
    max_abs_errs.reserve(samples);

    llps::utilities::constexpr_for<samples>([&]<size_t J>(llps::utilities::size_t_constant<J>)
    {
        static constexpr size_t rows = 16 + J * 10;
        static constexpr value_type dx = (x_max - x_min) / rows;

        using grid_t = llps::grid<value_type, rows, rows>;

        grid_t phi, expected;
        llps::apply_equi2D(phi, x_min, x_max, test_func::phi);
        llps::apply_equi2D(expected, x_min, x_max, test_func::dphi);

        grid_t actual = llps::calculus::laplacian_central_fd<TestFixture::error_order>(phi, dx, dx);
        //Using max absolute error to measure error
        value_type max_abs_err = llps::utilities::max_abs_error(expected, actual);
        
        delta_xs.push_back(dx);
        max_abs_errs.push_back(max_abs_err);
    });

    //Calculate machine imprecision point
    size_t stop_index = 1;
    while (stop_index < samples) {
        if (max_abs_errs[stop_index] >= max_abs_errs[stop_index - 1])
            break;
        else
            ++stop_index;
    }

    auto fit = llps::utilities::poly_fit1D(
        delta_xs     | std::ranges::views::take(stop_index),
        max_abs_errs | std::ranges::views::take(stop_index),
        static_cast<double (*)(double)>(std::log),
        static_cast<double (*)(double)>(std::log));

    ASSERT_GE(fit.gradient, expected);
}