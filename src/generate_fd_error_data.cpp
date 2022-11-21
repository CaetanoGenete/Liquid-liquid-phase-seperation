#include "grid.hpp"

#include <iostream>
#include <cmath>
#include <numbers>

int main()
{
    grid<double, 8, 8> test;
    apply_equi2D(test, 0., 2.*std::numbers::pi, [](double x, double y) {
        return std::cos(x) + std::sin(y);
    });

    for (size_t row = 0; row < test.rows(); ++row) {
        for (size_t col = 0; col < test.cols(); ++col)
            std::cout << test(row, col) << ", ";
        std::cout << std::endl;
    }
}