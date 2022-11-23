#ifndef GRID_HPP_INCLUDED
#define GRID_HPP_INCLUDED

#include <vector>      //Access to std::vector
#include <type_traits> //Access to std::is_same_v

#include "aligned_allocator.hpp"

#ifdef LLPS_USE_MKL
    template<class Type>
    using _grid_default_alloc = fftw_allocator<Type>;
#else
    template<class Type>
    using _grid_default_alloc = std::allocator<Type>;
#endif // MKL_USE_MKL


template<typename Type, size_t _rows, size_t _cols, typename Container = std::vector<Type, _grid_default_alloc<Type>>>
struct grid;

template<typename Type, size_t _rows, size_t _cols, typename Allocator>
struct grid<Type, _rows, _cols, std::vector<Type, Allocator>>
{
private:
    using _underlying_t = std::vector<Type, Allocator>;

public:
    using value_type      = typename _underlying_t::value_type;
    using reference       = typename _underlying_t::reference;
    using const_reference = typename _underlying_t::const_reference;
    using size_type       = typename _underlying_t::size_type;

    using iterator        = typename _underlying_t::iterator;
    using const_iterator  = typename _underlying_t::const_iterator;

    static_assert(std::is_same_v<Type, value_type>, "Container type mismatch!");

public:
    constexpr grid(): _underlying(size()) {}

public:
    constexpr const_reference operator()(size_type row, size_type column) const
    {
        //Row major order
        return _underlying[column + row * _cols];
    }

    constexpr reference operator()(size_type row, size_type column)
    {
        return const_cast<reference>(static_cast<const grid&>(*this)(row, column));
    }

public:
    constexpr iterator begin()              { return _underlying.begin(); };
    constexpr const_iterator begin()  const { return _underlying.cbegin(); };
    constexpr const_iterator cbegin() const { return _underlying.begin(); };

    constexpr iterator end()              { return _underlying.end(); };
    constexpr const_iterator end()  const { return _underlying.cend(); };
    constexpr const_iterator cend() const { return _underlying.end(); };

public:
    static constexpr size_type size() noexcept { return _rows * _cols; }
    static constexpr size_type rows() noexcept { return _rows; }
    static constexpr size_type cols() noexcept { return _cols; }

public:
    constexpr       value_type* data()       { return _underlying.data(); }
    constexpr const value_type* data() const { return _underlying.data(); }

private:
    _underlying_t _underlying;
};

template<typename Type, size_t _rows, size_t _cols, typename Allocator, class Callable>
constexpr void apply_equi2D(
    grid<Type, _rows, _cols, Allocator>& grid, 
    const Type x_min, const Type x_max, 
    const Type y_min, const Type y_max, 
    Callable func)
{
    const Type dx = (x_max - x_min) / _cols;
    const Type dy = (y_max - y_min) / _rows;

    for (size_t row = 0; row < _rows; ++row)
        for (size_t col = 0; col < _cols; ++col)
            grid(row, col) = std::invoke(func, col * dx, row * dy);
}

template<typename Type, size_t _rows, size_t _cols, typename Allocator, class Callable>
constexpr void apply_equi2D(
    grid<Type, _rows, _cols, Allocator>& grid,
    const Type x_min, const Type x_max,
    Callable func)
{
    return apply_equi2D(grid, x_min, x_max, x_min, x_max, func);
}


#endif // !GRID_HPP_INCLUDED