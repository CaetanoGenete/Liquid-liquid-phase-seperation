#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>

#include "boost/numeric/odeint.hpp"

#include "_modelb_common.hpp"

#include "llps/utilities/io.hpp"
#include "llps/utilities/timer.hpp"
#include "llps/calculus/differentiate.hpp"
#include "llps/grid.hpp"


template<size_t rows, size_t cols>
class coupled_field : public llps::grid<double, rows*2, cols>
{
private:
    using _base_t = llps::grid<double, rows * 2, cols>;

public:
    using _base_t::_base_t;

public:
    auto field(size_t index)
    {
        return llps::subgrid_view<_base_t, rows, cols>(*this, index * rows, 0);
    }

    auto field(size_t index) const
    {
        return llps::const_subgrid_view<_base_t, rows, cols>(*this, index * rows, 0);
    }
};

using state_type = coupled_field<256, 256>;

template<size_t order>
struct modelb_coupled
{
public:
    modelb_coupled(double a, double b, double k) :
        _a(a), _b(b), _k(k), _temp() {}

public:
    void operator()(const state_type& phi, state_type& dphi, double)
    {
        
        static constexpr double dx = 1.;
        static constexpr double dy = 1.;

        static constexpr double xi[] = {2., -1.};

        static constexpr size_t field_size = state_type::size() / 2;

        /*
        llps::calculus::laplacian_central_fd<order>(phi, dphi, dx, dy);

        auto dphi_it = dphi.begin();
        for (auto phi_it = phi.begin(); phi_it != phi.end(); ++phi_it, ++dphi_it) {
            const auto& phi = *phi_it;
            *dphi_it = phi * (_a + _b * phi * phi) - _k * (*dphi_it);
        }

        dphi = llps::calculus::laplacian_central_fd<order>(dphi, dx, dy);
        */

        for (size_t i = 0; i < 2; ++i) {
            auto temp_field = _temp.field(i);
            llps::calculus::laplacian_central_fd<order>(phi.field(i), temp_field, dx, dy);
        }

        auto temp_it = _temp.begin();
        auto phi_it = phi.begin();
        auto other_it = phi_it + field_size;
        for (size_t j = 0; j < 2; ++j) {
            for (size_t i = 0; i < field_size; ++i, ++temp_it, ++phi_it, ++other_it)
            {
                const auto& phi = *phi_it;
                *temp_it = phi * (_a + _b * phi * phi) - _k * (*temp_it) + xi[j] * (*other_it);
            }

            other_it = phi.begin();
        }

        /*

        for (auto phi_it = phi.begin(); phi_it != phi.end(); ++phi_it, ++temp_it) {
            const auto& phi = *phi_it;
            *temp_it = phi * (_a + _b * phi * phi) - _k * (*temp_it);
        }
        */

        for (size_t i = 0; i < 2; ++i) {
            auto dphi_field = dphi.field(i);
            llps::calculus::laplacian_central_fd<order>(_temp.field(i), dphi_field, dx, dy);
        }
    }

private:
    double _a, _b, _k;
    state_type _temp;
};


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
    constexpr double dt = 1.;

    //Sampling 
    constexpr double sample_int = 1.;
    constexpr size_t samples = static_cast<size_t>((t_max - t_min) / sample_int) + 1;

    using field_type = llps::grid<state_type::value_type, state_type::rows() / 2, state_type::cols()>;
    std::vector<field_type> field1;
    std::vector<field_type> field2;
    field1.reserve(samples);
    field2.reserve(samples);

    auto model = modelb_coupled<6>(a, b, k);

    { llps::timer timer;

    double last_t = t_min;
    odeint::integrate_adaptive(stepper, std::ref(model), phi0, t_min, t_max, dt, [&](const state_type& phi, double t) {
        if (t - last_t >= sample_int) {
            std::cout << "Progress: " << std::setprecision(2) << t << "/" << std::fixed << t_max << "\r";

            field1.push_back(phi.field(0));
            field2.push_back(phi.field(1));
            last_t += sample_int;
        }
        });
    }

    save_to_file(LLPS_OUTPUT_DIR"modelb_coupled_1(a=-b=-k=-1).dat", field1, "$\\phi_1$");
    save_to_file(LLPS_OUTPUT_DIR"modelb_coupled_2(a=-b=-k=-1).dat", field2, "$\\phi_2$");

}