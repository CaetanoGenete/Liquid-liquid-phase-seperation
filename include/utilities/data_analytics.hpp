#ifndef UTILS_DATA_ANALYTICS_HPP_INCLUDED
#define UTILS_DATA_ANALYTICS_HPP_INCLUDED

#include <iterator>
#include <assert.h>

namespace utilities {

    template<typename Type>
    struct linear_regression_pair
    {
        Type intercept;
        Type gradient;
    };

    template<typename XIt, typename YIt>
    inline constexpr auto poly_fit1D(XIt first_x, XIt last_x, YIt first_y, YIt last_y)
    {
        using x_value = std::iter_value_t<XIt>;
        using y_value = std::iter_value_t<YIt>;

        size_t size = 0;

        x_value x_bar{};
        y_value y_bar{};

        XIt x_it = first_x;
        YIt y_it = first_y;
        for (; x_it != last_x && y_it != last_y; ++x_it, ++y_it, ++size) {
            x_bar += *x_it;
            y_bar += *y_it;
        }

        //Calculate means
        x_bar /= size;
        y_bar /= size;

        //Ranges must be of the same length
        assert(x_it == last_x && y_it == last_y);

        x_it = first_x;
        y_it = first_y;

        y_value num{};
        x_value den{};
        for (; x_it != last_x; ++x_it, ++y_it)
        {
            x_value x_diff = (*x_it - x_bar);

            num += x_diff * (*y_it - y_bar);
            den += x_diff * x_diff;
        }

        y_value gradient = num / den;
        return linear_regression_pair{ y_bar - gradient * x_bar, gradient };
    }
}

#endif // !UTILS_DATA_ANALYTICS_HPP_INCLUDED