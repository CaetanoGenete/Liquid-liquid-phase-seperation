
#include <iostream>
#include <array>
#include <random>
#include <algorithm>
#include <ranges>
#include <fstream>

#include "boost/numeric/odeint.hpp"

//#include "_modelb_common.hpp"
#include "multi_range_algebra.hpp"

#include "llps/grid.hpp"
#include "llps/utilities/timer.hpp"
#include "llps/utilities/io.hpp"
#include "llps/calculus/differentiate.hpp"

using time_type = double;

template<size_t order, typename _field_type>
struct modelb_coupled
{
private:
    using _state_type = std::array<_field_type, 2>;
    using _value_type = typename _field_type::value_type;

public:
    modelb_coupled(double a, double b, double k, std::array<double, 2> xi) :
        _a(a), _b(b), _k(k), _xi(xi) {}

public:
    LLPS_FORCE_INLINE void operator()(const _state_type& phi, _state_type& dphi, time_type)
    {
        static constexpr _value_type dx = 1.;
        static constexpr _value_type dy = 1.;

        for (size_t i = 0, j = 1; i < 2; j = i++) {
            const auto& field_i = phi[i];
            const auto& field_j = phi[j];
            auto& dfield_i = dphi[i];

            llps::calculus::laplacian_central_fd<order>(field_i, dfield_i, dx, dy);

            auto field_i_it = field_i.begin();
            auto field_j_it = field_j.begin();
            for (auto& dfield_i_val : dfield_i) {
                const _value_type field_i_val = *field_i_it;

                dfield_i_val = field_i_val * (_a + _b * field_i_val * field_i_val) - _k * dfield_i_val + _xi[i] * (*field_j_it);

                ++field_i_it;
                ++field_j_it;
            }

            dfield_i = llps::calculus::laplacian_central_fd<order>(dfield_i, dx, dy);
        }
    }

private:
    _value_type _a, _b, _k;
    std::array<_value_type, 2> _xi;
};

int main()
{
    using namespace boost::numeric;

    using field_type = llps::grid<double, 256, 256>;
    using state_type = std::array<field_type, 2>;
    using value_type = field_type::value_type;

    using stepper_type = odeint::runge_kutta_cash_karp54<state_type, value_type, state_type, time_type, array_of_ranges_algebra>;
    auto stepper = odeint::make_controlled<stepper_type>(1e-10, 1e-6);

    std::default_random_engine rnd_eng{ 69 };
    std::normal_distribution normal_dist{ 0., 1. };

    state_type phi0;
    value_type intphi0 = 0;
    for (auto& field : phi0) {
        std::ranges::generate(field, std::bind(normal_dist, rnd_eng));
        
        intphi0 = std::accumulate(field.begin(), field.end(), intphi0);
    }

    //Model B paramaters
    constexpr value_type a = -1.;
    constexpr value_type b = -a;
    constexpr value_type k = 1.;
    constexpr value_type xi1 = 2.;
    constexpr value_type xi2 = -1.;

    //Integration paramaters
    constexpr time_type t_min = 0.;
    constexpr time_type t_max = 1000.;
    constexpr time_type dt = 1.;

    //Sampling  parameters
    constexpr size_t samples = 1001;
    constexpr time_type sample_int = (t_max - t_min) / (samples - 1);

    std::vector<time_type> times;
    std::vector<field_type> fields[2];

    //Initialise outputs
    times.reserve(samples);

    for(std::vector<field_type>& field : fields)
        field.reserve(samples);

    //Integration:

    auto model = modelb_coupled<6, field_type>(a, b, k, {xi1, xi2});
    { 
        llps::timer timer;

        time_type last_t = -sample_int;
        odeint::integrate_adaptive(stepper, model, phi0, t_min, t_max, dt, [&](const state_type& phi, time_type t)
        {
            if (t - last_t >= sample_int) {
                std::cout << "Progress: " << std::setprecision(2) << t << "/" << std::fixed << t_max << "\r";

                for (size_t i = 0; i < std::ranges::size(fields); ++i)
                    fields[i].push_back(phi[i]);

                times.push_back(t);
                last_t += sample_int;
            }
        });
    }

    value_type vmin = std::numeric_limits<value_type>::max();
    value_type vmax = std::numeric_limits<value_type>::min();
    for (const auto& field : fields) {
        for (const auto& frame : field) {
            auto [frame_min, frame_max] = std::ranges::minmax(frame);
            vmin = std::min(vmin, frame_min);
            vmax = std::max(vmax, frame_max);
        }
    }

    std::string data_suffix = std::format(
        "phi0={:.2f},a={:.2f},b={:.2f},k={:.2f},xi_1={:.2f},xi_2={:.2f},t={}", intphi0, a, b, k, xi1, xi2, t_max);

    std::ofstream file(LLPS_OUTPUT_DIR"simulations/coupled model B/coupled_modelB(" + data_suffix + ").dat", std::ios::binary);

    std::string title_data = std::format(
        "$\\phi_0$={:.2f}, a={:.2f}, b={:.2f}, $\\kappa$={:.2f}, $\\xi_1$={:.2f}, $\\xi_2$={:.2f}", intphi0, a, b, k, xi1, xi2);

    llps::utilities::plot_header plot_header;
    plot_header.title = "Coupled ModelB fields ($\\phi_1, \\phi_2$). With params:\n" + title_data;
    plot_header.x_label = "x";
    plot_header.y_label = "y";

    llps::utilities::serialise_plot_header(file, 2, plot_header);

    llps::utilities::serialise_meta_data_header(file, 2);
    llps::utilities::serialise_meta_data(file, "vmin", vmin);
    llps::utilities::serialise_meta_data(file, "vmax", vmax);

    for (size_t i = 0; i < 2; ++i) {
        llps::utilities::video_header<value_type, value_type> video_header;
        video_header.sub_title = "$\\phi_" + std::to_string(i+1) + "$";

        llps::utilities::serialise_video_header(file, field_type::cols(), field_type::rows(), fields[i].size(), video_header);

        for (auto& frame : fields[i]) {
            for (auto& value : frame)
                llps::utilities::serialise_to_binary(file, value);
        }
    }

    for (auto& time : times)
        llps::utilities::serialise_to_binary(file, time);

    file.close();

    //Pause
    std::cout << "Saving succeeded!";
    std::cin.get();
}