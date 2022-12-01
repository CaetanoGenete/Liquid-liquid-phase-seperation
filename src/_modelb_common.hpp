#ifndef _MODELB_COMMON_HPP_INCLUDED
#define _MODELB_COMMON_HPP_INCLUDED

#include <vector>
#include <algorithm>
#include <fstream>
#include <string>

#include "calculus/differentiate.hpp"
#include "utilities/io.hpp"
#include "grid.hpp"

using state_type = llps::grid<double, 256, 256>;

template<size_t order>
struct modelb
{
public:
    modelb(double a, double b, double k) :
        _a(a), _b(b), _k(k) {}

public:
    void operator()(const state_type& phi, state_type& dphi, double)
    {
        static constexpr double dx = 1.;
        static constexpr double dy = 1.;
        llps::calculus::laplacian_central_fd<order>(phi, dphi, dx, dy);

        auto dphi_it = dphi.begin();
        for (auto phi_it = phi.begin(); phi_it != phi.end(); ++phi_it, ++dphi_it) {
            const auto& phi = *phi_it;
            *dphi_it = phi * (_a + _b * phi * phi) - _k * (*dphi_it);
        }

        dphi = llps::calculus::laplacian_central_fd<order>(dphi, dx, dy);
    }

private:
    double _a, _b, _k;
};

void save_to_file(const char* file_name, const std::vector<state_type>& data, std::string title)
{
    //Calculate minimum and maximum values
    auto [vmin, vmax] = std::ranges::minmax(data.front());
    for (auto& frame : data) {
        auto [phi_min, phi_max] = std::ranges::minmax(frame);
        vmin = std::min(phi_min, vmin);
        vmax = std::max(phi_max, vmax);
    }

    std::ofstream file(file_name, std::ios::binary);

    llps::utilities::plot_header plot_header;
    plot_header.title = title;
    plot_header.x_label = "x";
    plot_header.y_label = "y";

    llps::utilities::serialise_plot_header(file, 1, plot_header);

    llps::utilities::video_header<state_type::value_type, state_type::value_type> video_header;
    video_header.min_max = { vmin, vmax };

    llps::utilities::serialise_video_header(file, state_type::cols(), state_type::rows(), data.size(), video_header);

    for (auto& frame : data) {
        for (auto& value : frame)
            llps::utilities::serialise_to_binary(file, value);
    }

    file.close();
}

#endif // !_MODELB_COMMON_HPP_INCLUDED