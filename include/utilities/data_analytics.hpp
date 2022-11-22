#ifndef UTILS_DATA_ANALYTICS_HPP_INCLUDED
#define UTILS_DATA_ANALYTICS_HPP_INCLUDED

#include <iterator> //Access to std::iter_value_t and iterator concepts
#include <assert.h> //Access to assert macro
#include <cmath>    //Access to std::abs

namespace utilities {

    template<typename Type>
    struct linear_regression_pair
    {
        Type intercept;
        Type gradient;
    };

    template<typename XIt, typename YIt>
    constexpr auto poly_fit1D(XIt first_x, XIt last_x, YIt first_y, YIt last_y)
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

    template<std::input_iterator InputIt1, std::input_iterator InputIt2>
    constexpr auto max_abs_error(InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2)
    {
        auto result = std::abs(*first1 - *first2);
        ++first1;
        ++first2;

        for (; first1 != last1 && first2 != last2; ++first1, ++first2)
            result = std::max(result, std::abs(*first1 - *first2));

        //Function is undefined for ranges of different size
        assert(first1 == last1 && first2 == last2);

        return result;
    }
}

#endif // !UTILS_DATA_ANALYTICS_HPP_INCLUDED