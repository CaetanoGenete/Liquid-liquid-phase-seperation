#ifndef UTILITIES_META_HPP_INCLUDED
#define UTILITIES_META_HPP_INCLUDED

#include <utility>     //For access to std::index_sequence
#include <functional>  //For access to std::invoke
#include <type_traits> //For access to std::integral_constant

namespace utilities {

    template<size_t I>
    using size_t_constant = std::integral_constant<size_t, I>;

    template<size_t from, size_t step, size_t ... indicies, typename Callable>
    void _constexpr_for_helper(Callable&& callable, std::index_sequence<indicies...>)
    {
        (std::invoke(callable, size_t_constant<indicies* step + from>{}), ...);
    }

    template<size_t from, size_t to, size_t step = 1, class Callable>
    void constexpr_for(Callable&& callable)
    {
        _constexpr_for_helper<from, step>(std::forward<Callable>(callable), std::make_index_sequence<(to - from + step - 1) / step>{});
    }

    template<size_t to, class Callable>
    void constexpr_for(Callable&& callable)
    {
        constexpr_for<0, to>(std::forward<Callable>(callable));
    }

}

#endif // !UTILITIES_META_HPP_INCLUDED