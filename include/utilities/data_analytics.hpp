#ifndef UTILS_DATA_ANALYTICS_HPP_INCLUDED
#define UTILS_DATA_ANALYTICS_HPP_INCLUDED

#include <iterator>

namespace utils {

    template<typename Type>
    struct linear_regression_pair
    {
        Type intercept;
        Type gradient;
    };
}

#endif // !UTILS_DATA_ANALYTICS_HPP_INCLUDED