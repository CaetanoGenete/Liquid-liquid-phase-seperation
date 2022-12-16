#include "llps/grid.hpp"

#include <iostream>
#include <array>
#include <random>
#include <algorithm>
#include <ranges>
#include <fstream>

#include "boost/numeric/odeint.hpp"
//#include "_modelb_common.hpp"
#include "multi_range_algebra.hpp"
#include "llps/utilities/timer.hpp"
#include "llps/utilities/io.hpp"
#include "llps/calculus/differentiate.hpp"

using time_type = double;


template<size_t order, typename _field_type>
struct modelb_coupled_diffusion
{
private:
    using _state_type = std::array<_field_type, 2>;
    using _value_type = typename _field_type::value_type;

public:
    modelb_coupled_diffusion(_value_type a, _value_type b, _value_type k, _value_type k01, _value_type k10, _value_type d):
        _a(a), _b(b), _k(k), _k01(k01), _k10(k10), _d(d) {}

public:
    LLPS_FORCE_INLINE void operator()(const _state_type& phi, _state_type& dphi, time_type)
    {
        static constexpr _value_type dx = 1.;
        static constexpr _value_type dy = 1.;

        const auto& field_1 = phi[0];
        auto& dfield_1 = dphi[0];

        llps::calculus::laplacian_central_fd<order>(field_1, dfield_1, dx, dy);

        auto field_1_it = field_1.begin();
        for (auto& dfield_1_val : dfield_1) {
            const _value_type field_1_val = *field_1_it;
            ++field_1_it;

            dfield_1_val = field_1_val * (_a + _b * field_1_val * field_1_val) - _k * dfield_1_val; //+ _xi[i] * (*field_j_it);
        }

        dfield_1 = llps::calculus::laplacian_central_fd<order>(dfield_1, dx, dy);

        //Diffusion
        llps::calculus::laplacian_central_fd<order>(phi[1], dphi[1], dx, dy);

             field_1_it = field_1.begin();
        auto field_2_it = phi[1].begin();

        auto dfield_1_it = dfield_1.begin();
        auto dfield_2_it = dphi[1].begin();

        for (; field_1_it != field_1.end(); ++field_1_it, ++field_2_it, ++dfield_1_it, ++dfield_2_it)
        {
            //Multiply diffusion coeff here cause why not. 
            *dfield_2_it *= _d;

            const _value_type switching = _k01 * (*field_2_it) -_k10 * (*field_1_it);
            *dfield_1_it += switching;
            *dfield_2_it -= switching;
        }
    }

private:
    _value_type _a, _b, _k, _k01, _k10, _d;
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
    std::normal_distribution normal_dist{ -0.3, 1. };

    state_type phi0;
    value_type intphi0 = 0;
    for (auto& field : phi0) {
        std::ranges::generate(field, std::bind(normal_dist, rnd_eng));

        intphi0 = std::accumulate(field.begin(), field.end(), intphi0);
    }
    //Ensure the fields are different
    assert(!std::ranges::equal(phi0[0], phi0[1]));

    //Model B paramaters
    constexpr value_type a = -1.;
    constexpr value_type b = -a;
    constexpr value_type k = 1.;
    constexpr value_type k01 = 0.;
    constexpr value_type k10 = 0.;

    //Diffusion parameters:
    constexpr value_type d = 1.;

    //Integration paramaters
    constexpr time_type t_min = 0.;
    constexpr time_type t_max = 1000;
    constexpr time_type dt = 1.;

    //Sampling  parameters
    constexpr size_t samples = 1001;
    constexpr time_type sample_int = (t_max - t_min) / (samples - 1);

    std::vector<time_type> times;
    std::vector<field_type> fields[2];

    //Initialise outputs
    times.reserve(samples);

    for (std::vector<field_type>& field : fields)
        field.reserve(samples);

    //Integration:
    auto model = modelb_coupled_diffusion<6, field_type>(a, b, k, k01, k10, d);
    {
        llps::timer timer;

        time_type last_t = 0;
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
        "phi0={:.2f},a={:.2f},b={:.2f},k={:.2f},k01={:.2E},k10={:.2E},D={:.2f},t={}", intphi0, a, b, k, k01, k10, d, t_max);

    std::ofstream file(LLPS_OUTPUT_DIR"simulations/coupled modelB/coupled_modelB_diffusion(" + data_suffix + ").dat", std::ios::binary);

    std::string title_data = std::format(
        "$\\phi_0$={:.2f}, a={:.2f}, b={:.2f}, $\\kappa$={:.2f}, $k_{{01}}$={:.2E}, $k_{{10}}$={:.2E}, D={:.2f}", intphi0, a, b, k, k01, k10, d);

    llps::utilities::plot_header plot_header;
    plot_header.title = "Coupled ModelB and diffusion field ($\\phi_1, \\phi_2$). With params:\n" + title_data;
    plot_header.x_label = "x";
    plot_header.y_label = "y";

    llps::utilities::serialise_plot_header(file, 2, plot_header);

    llps::utilities::serialise_meta_data_header(file, 2);
    llps::utilities::serialise_meta_data(file, "vmin", vmin);
    llps::utilities::serialise_meta_data(file, "vmax", vmax);

    for (size_t i = 0; i < 2; ++i) {
        llps::utilities::video_header<value_type, value_type> video_header;
        video_header.sub_title = "$\\phi_" + std::to_string(i + 1) + "$";

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