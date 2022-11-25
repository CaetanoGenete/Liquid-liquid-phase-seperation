#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <numbers>
#include <iostream>
#include <iomanip>

#include "grid.hpp"
#include "calculus/differentiate.hpp"
#include "utilities/data_analytics.hpp"
#include "utilities/meta.hpp"
#include "utilities/io.hpp"

//For access to s suffix
using namespace std::literals::string_literals;

struct test_func
{
    static double phi(double x, double y)
    {
        return std::exp(std::cos(x) + std::sin(y));
    }

    static double dphi(double x, double y)
    {

        return phi(x, y) * (std::cos(y) * std::cos(y) + std::sin(x) * std::sin(x) - std::sin(y) - std::cos(x));
    }
};

int main()
{
    //For pretty plots
    static constexpr const char* colours[] = {"#f0f921", "#fdb42f", "#ed7953", "#cc4778", "#9c179e", "#5c01a6", "#0d0887"};

    static constexpr double x_min = 0;
    static constexpr double x_max = 2.*std::numbers::pi;

    std::ofstream file(LLPS_OUTPUT_DIR"finite_difference_accuracy.dat", std::ios::binary);
 
    llps::utilities::plot_header plot_header;
    plot_header.title   = "log-log plot of max absolute error of $\\nabla^2 \\phi$ using finite differences.\n$\\phi(x, y) = e^{\\cos(x) + \\sin(x)}$";
    plot_header.x_label = "$\\Delta x = \\Delta y$";
    plot_header.y_label = "max absolute error";
    plot_header.x_scale = "log";
    plot_header.y_scale = "log";

    constexpr size_t lines_count = 7;
    llps::utilities::serialise_plot_header(file, lines_count, plot_header);

    llps::utilities::constexpr_for<lines_count>([&]<size_t I>(llps::utilities::size_t_constant<I>)
    {
        static constexpr size_t order = 2 * (I+1);
        static constexpr size_t samples = 30;

        std::vector<double> delta_xs;
        std::vector<double> max_abs_errs;
        delta_xs.reserve(samples);
        max_abs_errs.reserve(samples);

        llps::utilities::constexpr_for<samples>([&]<size_t J>(llps::utilities::size_t_constant<J>)
        {
            static constexpr size_t rows = 16 + J * 10;
            static constexpr double dx = (x_max - x_min) / rows;

            llps::grid<double, rows, rows> expected;
            llps::grid<double, rows, rows> phi;
            llps::apply_equi2D(expected, x_min, x_max, test_func::dphi);
            llps::apply_equi2D(phi, x_min, x_max, test_func::phi);

            auto actual = llps::calculus::laplacian_central_fd<order>(phi, dx, dx);
            
            double max_abs_err = llps::utilities::max_abs_error(expected.begin(), expected.end(), actual.begin(), actual.end());

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
            delta_xs.begin(), 
            delta_xs.begin() + stop_index,
            max_abs_errs.begin(), 
            max_abs_errs.begin() + stop_index,
            std::logl, std::logl);

        std::cout << "Actual error order: (expected: " << std::setw(2) << order << ") x^" << fit.gradient << std::endl;

        llps::utilities::line_header line_header;
        line_header.colour = colours[I];
        line_header.label  = "O($\\Delta x^{"s + std::to_string(order) + "}$)"s;

        llps::utilities::serialise_line_header<double, double>(file, samples, line_header);

        for (double dx : delta_xs)
            llps::utilities::serialise_to_binary(file, dx);
        for (double errors : max_abs_errs)
            llps::utilities::serialise_to_binary(file, errors);
    });

    file.close();
}