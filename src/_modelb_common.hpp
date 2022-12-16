#ifndef _MODELB_COMMON_HPP_INCLUDED
#define _MODELB_COMMON_HPP_INCLUDED

#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include <numeric>

#include "calculus/differentiate.hpp"
#include "utilities/io.hpp"
#include "grid.hpp"

//using state_type = llps::grid<double, 256, 256>;

template<size_t order, class state_type>
struct modelb
{
public:
    modelb(double a, double b, double k) :
        _a(a), _b(b), _k(k) {}

public:
    LLPS_FORCE_INLINE void operator()(const state_type& phi, state_type& dphi, double)
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

template<class FrameType>
void save_to_file(const char* file_name, const std::vector<FrameType>& data, std::string title)
{
    using value_type = FrameType::value_type;

    //Calculate minimum and maximum values
    auto vmin = std::numeric_limits<value_type>::max();
    auto vmax = std::numeric_limits<value_type>::min();
    for (auto& frame : data) {
        auto [frame_min, frame_max] = std::ranges::minmax(frame);
        vmin = std::min(frame_min, vmin);
        vmax = std::max(frame_max, vmax);
    }

    std::ofstream file(file_name, std::ios::binary);

    llps::utilities::plot_header plot_header;
    plot_header.title = title;
    plot_header.x_label = "x";
    plot_header.y_label = "y";

    llps::utilities::serialise_plot_header(file, 1, plot_header);

    llps::utilities::video_header<value_type, value_type> video_header;
    video_header.min_max = { vmin, vmax };

    llps::utilities::serialise_video_header(file, FrameType::cols(), FrameType::rows(), data.size(), video_header);

    for (auto& frame : data) {
        for (auto& value : frame)
            llps::utilities::serialise_to_binary(file, value);
    }

    file.close();
}

#endif // !_MODELB_COMMON_HPP_INCLUDED