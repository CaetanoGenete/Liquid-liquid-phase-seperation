#include "gtest/gtest.h"

#include <array>   //Access to std::array
#include <vector>  //Access to std::vector
#include <numeric> //Access to std::iota

#include "utilities/data_analytics.hpp"

TEST(analytics_tests, test_poly_fit1D)
{
    //Data acquired from: https://en.wikipedia.org/wiki/Simple_linear_regression

    constexpr auto x = std::array{ 1.47, 1.50, 1.52, 1.55, 1.57, 1.60, 1.63 ,1.65, 1.68, 1.70, 1.73, 1.75, 1.78, 1.80, 1.83 };
    constexpr auto y = std::array{ 52.21, 53.12, 54.48, 55.84, 57.20, 58.57, 59.93, 61.29, 63.11, 64.47, 66.28, 68.10, 69.92, 72.19, 74.46 };

    auto fit = utilities::poly_fit1D(x.begin(), x.end(), y.begin(), y.end());
    ASSERT_NEAR(fit.intercept, -39.062, 1e-4);
    ASSERT_NEAR(fit.gradient, 61.2722, 1e-4);
}

//Trival test case: x = y
TEST(analytics_tests, test_poly_fit1D_trivial)
{
    static constexpr size_t samples = 100;
    
    std::vector<double> x(samples);
    std::iota(x.begin(), x.end(), 0.);

    std::vector<double> y(x);

    auto fit = utilities::poly_fit1D(x.begin(), x.end(), y.begin(), y.end());
    ASSERT_EQ(fit.intercept, 0);
    ASSERT_EQ(fit.gradient, 1);
}
