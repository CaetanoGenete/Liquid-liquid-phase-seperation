#include "grid.hpp"

#include <iostream>

int main()
{
    using grid_t = llps::grid<double, 8, 8>;
    grid_t grid;
    llps::apply_equi2D(grid, 0, 8, [](double x, double y) { return x; });

    llps::subgrid_view<grid_t, 3, 3> test(grid, 2, 3);

    for (size_t row = 0; row < test.rows(); ++row) {
        for (size_t col = 0; col < test.cols(); ++col)
            test(row, col) = 2.;

    }

    for (size_t row = 0; row < grid.rows(); ++row) {
        for (size_t col = 0; col < grid.cols(); ++col)
            std::cout << grid(row, col) << ", ";

        std::cout << std::endl;
    }


}