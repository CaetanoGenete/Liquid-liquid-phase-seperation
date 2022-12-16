#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <numbers>
#include <iostream>
#include <iomanip>
#include <ranges>

#include "llps/grid.hpp"
#include "llps/calculus/differentiate.hpp"
#include "llps/utilities/data_analytics.hpp"
#include "llps/utilities/meta.hpp"
#include "llps/utilities/io.hpp"

//For access to s suffix
using namespace std::literals::string_literals;
//Use double precision for these plots
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

int main()
{
    //For pretty plots
    static constexpr const char* colours[] = {"#f0f921", "#fdb42f", "#ed7953", "#cc4778", "#9c179e", "#5c01a6", "#0d0887"};

    static constexpr value_type x_min = 0;
    static constexpr value_type x_max = 2.*std::numbers::pi;

    std::ofstream file(LLPS_OUTPUT_DIR"finite_difference_accuracy.dat", std::ios::binary);
 
    llps::utilities::plot_header plot_header;
    plot_header.title   = "log-log plot of max absolute error of $\\nabla^2 \\phi$ using finite differences.\n$\\phi(x, y) = e^{\\cos(x) + \\sin(x)}$";
    plot_header.x_label = "$\\Delta x = \\Delta y$";
    plot_header.y_label = "max absolute error";
    plot_header.x_scale = "log";
    plot_header.y_scale = "log";

    static constexpr size_t lines_count = 7;
    llps::utilities::serialise_plot_header(file, lines_count, plot_header);

    llps::utilities::constexpr_for<lines_count>([&]<size_t I>(llps::utilities::size_t_constant<I>)
    {
        static constexpr size_t order = 2 * (I+1);
        static constexpr size_t samples = 30;

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

            grid_t actual = llps::calculus::laplacian_central_fd<order>(phi, dx, dx);
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
            std::logl, std::logl);

        std::cout << "Actual error order: (expected: " << std::setw(2) << order << ") x^" << fit.gradient << std::endl;

        llps::utilities::line_header line_header;
        line_header.colour = colours[I];
        line_header.label  = "O($\\Delta x^{"s + std::to_string(order) + "}$)"s;

        llps::utilities::serialise_line_header<value_type, value_type>(file, samples, line_header);

        for (value_type dx : delta_xs)
            llps::utilities::serialise_to_binary(file, dx);
        for (value_type errors : max_abs_errs)
            llps::utilities::serialise_to_binary(file, errors);
    });

    file.close();
}