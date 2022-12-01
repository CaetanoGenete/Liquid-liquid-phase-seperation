#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>

#include "boost/numeric/odeint.hpp"

#include "_modelb_common.hpp"

#include "utilities/io.hpp"
#include "utilities/timer.hpp"
#include "calculus/differentiate.hpp"
#include "grid.hpp"

using state_type = llps::grid<double, 256, 256>;

int main()
{
    using namespace boost::numeric;

    using stepper_type = odeint::runge_kutta_cash_karp54<state_type, state_type::value_type>;
    auto stepper = odeint::make_controlled<stepper_type>(1e-10, 1e-6);

    std::default_random_engine rnd_eng{ 69 };
    std::normal_distribution normal_dist{ 0., 1. };

    state_type phi0;
    std::ranges::generate(phi0, std::bind(normal_dist, rnd_eng));

    //Model B paramaters
    constexpr double a = -1.;
    constexpr double b = -a;
    constexpr double k = 1.;

    //Integration paramaters
    constexpr double t_min = 0.;
    constexpr double t_max = 1000.;

    constexpr double sample_int = 1.;
    constexpr size_t samples = static_cast<size_t>((t_max - t_min)/sample_int) + 1;
    constexpr double dt = 1.;

    std::vector<state_type> result;
    result.reserve(samples);

    { llps::timer timer;

        double last_t = t_min;
        odeint::integrate_adaptive(stepper, modelb<6>(a, b, k), phi0, t_min, t_max, dt, [&](const state_type& phi, double t) {
            if (t - last_t >= sample_int) {
                std::cout << "Progress: " << std::setprecision(2) << t << "/" << std::fixed << t_max << "\r";
                result.push_back(phi);
                last_t += sample_int;
            }
        });
    }

    save_to_file(LLPS_OUTPUT_DIR"modelb(a=-b=-k=-1).dat", result, "Modelb simulation using finite difference,\nup to t=" + std::to_string(t_max));
}