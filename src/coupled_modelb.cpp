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

struct coupled_field:
    boost::additive1< coupled_field,
    boost::additive2< coupled_field, double,
    boost::multiplicative2< coupled_field, double > > >
{
public:
    coupled_field(): _fields() {}

    coupled_field& operator+=(const coupled_field& other)
    {
        for (char i = 0; i < std::ranges::size(_fields); ++i)
        {
            auto& other_field = other._fields[i];

            auto this_it = _fields[i].begin();
            for (auto other_it = other_field.begin(); other_it != other_field.end(); ++other_it, ++this_it)
                *this_it += *other_it;
        }

        return *this;
    }

    coupled_field& operator+=(const double value)
    {
        for (auto& field : _fields)
        {
            for (auto it = field.begin(); it != field.end(); ++it)
                *it += value;
        }

        return *this;
    }

    coupled_field& operator*=(const double scalar)
    {
        for (auto& field : _fields)
        {
            for (auto it = field.begin(); it != field.end(); ++it)
                *it *= scalar;
        }

        return *this;
    }

    coupled_field& operator/=(const coupled_field& other)
    {
        for (char i = 0; i < std::ranges::size(_fields); ++i)
        {
            auto& other_field = other._fields[i];

            auto this_it = _fields[i].begin();
            for (auto other_it = other_field.begin(); other_it != other_field.end(); ++other_it, ++this_it)
                *this_it /= *other_it;
        }

        return *this;
    }

    friend coupled_field operator/(const coupled_field& p1, const coupled_field& p2)
    {
        coupled_field result;

        for (char i = 0; i < std::ranges::size(result._fields); ++i)
        {
            auto& p1_field = p1._fields[i];
            auto& p2_field = p2._fields[i];

            auto result_it = result._fields[i].begin();
            auto p1_it = p1_field.begin();
            auto p2_it = p2_field.begin();
            for (; p1_it != p1_field.end(); ++p1_it, ++p2_it, ++result_it)
                *result_it = *p1_it / *p2_it;
        }

        return result;
    }

    friend coupled_field abs(const coupled_field& fields)
    {
        coupled_field result;

        for (char i = 0; i < std::ranges::size(result._fields); ++i)
        {
            auto& field = fields._fields[i];

            auto result_it = result._fields[i].begin();
            for (auto field_it = field.begin(); field_it != field.end(); ++field_it, ++result_it)
                *result_it = std::abs(*field_it);
        }

        return result;
    }

public:
    const state_type& field(char index) const noexcept { return _fields[index]; }

public:
    state_type _fields[2];
};

namespace boost::numeric::odeint {
    template<>
    struct vector_space_norm_inf< coupled_field >
    {
        typedef double result_type;
        double operator()(const coupled_field& fields) const
        {   
            return std::max(
                std::ranges::max(fields.field(0), {}, std::abs<double>),
                std::ranges::max(fields.field(1), {}, std::abs<double>));
        }
    };
}


template<size_t order>
struct coupled_modelb
{
public:
    coupled_modelb(double a, double b, double k) :
        _a(a), _b(b), _k(k) {}

public:
    void operator()(const coupled_field& phi, coupled_field& dphi, double)
    {
        static constexpr double dx = 1.;
        static constexpr double dy = 1.;
        
        for (char i = 0; i < std::ranges::size(phi._fields); ++i) {
            auto& phi_field = phi._fields[i];
            auto& dphi_field = dphi._fields[i];

            llps::calculus::laplacian_central_fd<order>(phi_field, dphi_field, dx, dy);

            auto dphi_it = dphi_field.begin();
            for (auto phi_it = phi_field.begin(); phi_it != phi_field.end(); ++phi_it, ++dphi_it) {
                const auto& phi = *phi_it;
                *dphi_it = phi * (_a + _b * phi * phi) - _k * (*dphi_it);
            }

            dphi_field = llps::calculus::laplacian_central_fd<order>(dphi_field, dx, dy);
        }
    }

private:
    double _a, _b, _k;
};

int main()
{
    using namespace boost::numeric;

    using stepper_type = odeint::runge_kutta_cash_karp54<coupled_field, double, coupled_field, double, odeint::vector_space_algebra>;
    auto stepper = odeint::make_controlled<stepper_type>(1e-10, 1e-6);

    std::default_random_engine rnd_eng{ 69 };
    std::normal_distribution normal_dist{ 0., 1. };

    coupled_field phi0;
    std::ranges::generate(phi0._fields[0], std::bind(normal_dist, rnd_eng));
    std::ranges::generate(phi0._fields[1], std::bind(normal_dist, rnd_eng));

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

    std::vector<coupled_field> result;
    result.reserve(samples);

    { llps::timer timer;

    double last_t = t_min;
    odeint::integrate_adaptive(stepper, coupled_modelb<6>(a, b, k), phi0, t_min, t_max, dt, [&](const coupled_field& phi, double t) {
        if (t - last_t >= sample_int) {
            std::cout << "Progress: " << std::setprecision(2) << t << "/" << std::fixed << t_max << "\r";

            result.push_back(phi);
            last_t += sample_int;
        }
        });
    }

    //save_to_file(LLPS_OUTPUT_DIR"modelb(a=-b=-k=-1).dat", result, "Modelb simulation using finite difference,\nup to t=" + std::to_string(t_max));
}