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
    using value_type = double;
    using fftw_vector_t = std::vector<value_type, llps::fftw_allocator<value_type>>;

    //For pretty plots
    static constexpr const char* colours[] = { "#f0f921", "#fdb42f", "#ed7953", "#cc4778", "#9c179e", "#5c01a6", "#0d0887" };

    static constexpr value_type x_min = 0;
    static constexpr value_type x_max = 2. * std::numbers::pi;

    std::ofstream file(LLPS_OUTPUT_DIR"spectral_accuracy.dat", std::ios::binary);

    llps::utilities::plot_header plot_header;
    plot_header.title   = "Plot of max absolute error of $\\nabla^2 \\phi$ using fourier spectral method.\n$\\phi(x, y) = e^{\\cos(x) + \\sin(x)}$";
    plot_header.x_label = "$\\Delta x = \\Delta y$";
    plot_header.y_label = "max absolute error";
    plot_header.x_scale = "log";
    plot_header.y_scale = "log";

    constexpr size_t lines_count = 1;
    llps::utilities::serialise_plot_header(file, lines_count, plot_header);

    static constexpr size_t samples = 10;

    std::vector<value_type> delta_xs;
    std::vector<value_type> max_abs_errs;
    delta_xs.reserve(samples);
    max_abs_errs.reserve(samples);

    llps::utilities::constexpr_for<samples>([&]<size_t I>(llps::utilities::size_t_constant<I>)
    {
        static constexpr size_t rows = 1 << (I + 2);
        static constexpr value_type dx = (x_max - x_min) / rows;

        llps::grid<value_type, rows, rows, fftw_vector_t> expected;
        llps::grid<value_type, rows, rows, fftw_vector_t> phi;
        llps::apply_equi2D(expected, x_min, x_max, test_func::dphi);
        llps::apply_equi2D(phi, x_min, x_max, test_func::phi);

        llps::calculus::laplacian_spectral(phi, dx, dx);

        value_type max_abs_err = llps::utilities::max_abs_error(expected.begin(), expected.end(), phi.begin(), phi.end());

        delta_xs.push_back(dx);
        max_abs_errs.push_back(max_abs_err);
    });

    llps::utilities::line_header line_header;
    line_header.colour = "red";

    llps::utilities::serialise_line_header<value_type, value_type>(file, samples, line_header);

    for (value_type dx : delta_xs)
        llps::utilities::serialise_to_binary(file, dx);
    for (value_type errors : max_abs_errs)
        llps::utilities::serialise_to_binary(file, errors);

    file.close();
}