#ifndef LLPS_GRID_HPP_INCLUDED
#define LLPS_GRID_HPP_INCLUDED

#include <vector>      //Access to std::vector
#include <type_traits> //Access to std::is_same_v

#include "aligned_allocator.hpp"

namespace llps {

#ifdef LLPS_USE_MKL
    template<class Type>
    using _grid_default_alloc = fftw_allocator<Type>;
#else
    template<class Type>
    using _grid_default_alloc = std::allocator<Type>;
#endif // LLPS_USE_MKL

    /*
    * Defines compile-time attributes of _basic_grid class. 
    */
    template<typename Type, size_t _rows, size_t _cols, typename Container>
    struct _grid_meta_data
    {
    public:
        using underlying_type = Container;

        using value_type      = typename underlying_type::value_type;
        using reference       = typename underlying_type::reference;
        using const_reference = typename underlying_type::const_reference;
        using size_type       = size_t;

        using iterator        = typename underlying_type::iterator;
        using const_iterator  = typename underlying_type::const_iterator;

        static_assert(std::is_same_v<Type, value_type>, "Container type mismatch!");

    public:
        constexpr static size_t rows = _rows;
        constexpr static size_t cols = _cols;
    };

    template<class Meta>
    struct _grid_base;

    template<typename Type, size_t _rows, size_t _cols, typename Container>
    struct _grid_base<_grid_meta_data<Type, _rows, _cols, Container>>
    {
    private:
        using _meta_type = _grid_meta_data<Type, _rows, _cols, Container>;

    public:
        using value_type      = typename _meta_type::value_type;
        using reference       = typename _meta_type::reference;
        using const_reference = typename _meta_type::const_reference;
        using size_type       = typename _meta_type::size_type;

        using iterator        = typename _meta_type::iterator;
        using const_iterator  = typename _meta_type::const_iterator;

    public:
        static consteval size_t size() noexcept { return _rows * _cols; }
        static consteval size_t rows() noexcept { return _rows; }
        static consteval size_t cols() noexcept { return _cols; }
    };

    /*
    * Two dimensional grid container with static size. Represents the underlying
    * container type as an mxn-grid whose points may be accessed through a row
    * index and a column index. 
    * 
    * Note: _basic_grid, by default, assumes Container is contiguous and satisfies
    * the following:
    * 
    * Container::value_type
    * Container::reference
    * Container::const_reference
    * Container::iterator
    * Container::const_iterator
    * 
    * Container(size_t)
    * Container::operator[](size_t) (const) -> Container::(const_)reference
    * Container::begin() (const) -> Container::(const_)iterator
    * Container::cbegin() const -> Container::const_iterator
    * Container::end() (const) -> Container::(const_)iterator
    * Container::cend() const -> Container::const_iterator
    * 
    * 
    * 
    */

    /*
    template<class Meta>
    struct _basic_grid;
    */

    template<class Meta>
    struct _basic_grid : public _grid_base<Meta>
    {
    private:
        using _base_t = _grid_base<Meta>;

    public:
        constexpr _basic_grid(): 
            _underlying(_base_t::size()) {}

    public:
        constexpr _base_t::const_reference operator()(_base_t::size_type row, _base_t::size_type column) const
        {
            //Row major order
            return _underlying[column + row * _base_t::cols()];
        }

        constexpr _base_t::reference operator()(_base_t::size_type row, _base_t::size_type column)
        {
            return const_cast<_base_t::reference>(static_cast<const _basic_grid&>(*this)(row, column));
        }

    public:
        constexpr _base_t::iterator begin()              { return _underlying.begin(); };
        constexpr _base_t::const_iterator begin()  const { return _underlying.cbegin(); };
        constexpr _base_t::const_iterator cbegin() const { return _underlying.begin(); };

        constexpr _base_t::iterator end()              { return _underlying.end(); };
        constexpr _base_t::const_iterator end()  const { return _underlying.cend(); };
        constexpr _base_t::const_iterator cend() const { return _underlying.end(); };

    public:
        constexpr       _base_t::value_type* data()       { return _underlying.data(); }
        constexpr const _base_t::value_type* data() const { return _underlying.data(); }

    private:
        typename Meta::underlying_type _underlying;
    };

    template<
        class Type,
        size_t _rows,
        size_t _cols,
        class Container = std::vector<Type, _grid_default_alloc<Type>>>
    using grid = _basic_grid<_grid_meta_data<Type, _rows, _cols, Container>>;

    template<class Type>
    struct grid_value;

    template<class Type>
    using grid_value_t = typename grid_value<Type>::type;

    template<typename Type, size_t _rows, size_t _cols, typename Container>
    struct grid_value<_grid_meta_data<Type, _rows, _cols, Container>>
    { using type = Type; };

    template<class Meta>
    struct grid_value<_basic_grid<Meta>>
    { using type = grid_value_t<Meta>; };


    template<typename Grid, size_t _view_rows, size_t _view_cols>
    struct const_subgrid_view;

    template<class GridMeta, size_t _view_rows, size_t _view_cols>
    struct const_subgrid_view<_basic_grid<GridMeta>, _view_rows, _view_cols>:
        public _grid_base<_grid_meta_data<grid_value_t<GridMeta>, _view_rows, _view_cols, typename GridMeta::underlying_type>>
    {
    public:
        using parent_grid_t = _basic_grid<GridMeta>;

    public:
        static_assert(_view_cols <= GridMeta::cols && _view_rows <= GridMeta::rows,
            "subgrid cannot exceed dimensions of its parent grid");

    public:
        constexpr const_subgrid_view(const parent_grid_t& grid, const_subgrid_view::size_type row_offset = 0, const_subgrid_view::size_type col_offset = 0) noexcept:
            _grid(grid),
            _row_offset(row_offset),
            _col_offset(col_offset)
        {}

        const_subgrid_view(const const_subgrid_view&) = delete;
        const_subgrid_view(const_subgrid_view&&) = delete;

    public:
        const_subgrid_view& operator=(const const_subgrid_view&) = delete;
        const_subgrid_view& operator=(const_subgrid_view&&) = delete;

    public:
        constexpr const_subgrid_view::const_reference operator()(const_subgrid_view::size_type row, const_subgrid_view::size_type col) const {
            return _grid(row + _row_offset, col + _col_offset); 
        }

    protected:
        const parent_grid_t& _grid;

        size_t _row_offset;
        size_t _col_offset;
    };

    template<typename Type, size_t _rows, size_t _cols, typename Container>
    const_subgrid_view(const grid<Type, _rows, _cols, Container>&) ->const_subgrid_view<grid<Type, _rows, _cols, Container>, _rows, _cols>;

    template<class Grid, size_t _view_rows, size_t _view_cols>
    struct subgrid_view :
        public const_subgrid_view<Grid, _view_rows, _view_cols>
    {
    private:
        using _base_t = const_subgrid_view<Grid, _view_rows, _view_cols>;

    public:
        constexpr subgrid_view(Grid& grid, _base_t::size_type row_offset = 0, _base_t::size_type col_offset = 0) noexcept :
            _base_t(grid, row_offset, col_offset) {}

    public:
        constexpr _base_t::reference operator()(_base_t::size_type row, _base_t::size_type col)
        { 
            //Cast is safe here since constructor requires non-const grid
            return const_cast<Grid&>(_base_t::_grid)(row + _base_t::_row_offset, col + _base_t::_col_offset);
        }
    };

    template<typename Type, size_t _rows, size_t _cols, typename Container>
    subgrid_view(grid<Type, _rows, _cols, Container>&)->subgrid_view<grid<Type, _rows, _cols, Container>, _rows, _cols>;


    template<typename Meta, typename Callable>
    constexpr void apply_equi2D(
        _basic_grid<Meta>& grid,
        const grid_value_t<Meta> x_min, const grid_value_t<Meta> x_max,
        const grid_value_t<Meta> y_min, const grid_value_t<Meta> y_max,
        Callable func)
    {
        using value_type = grid_value_t<Meta>;

        const value_type dx = (x_max - x_min) / grid.cols();
        const value_type dy = (y_max - y_min) / grid.rows();

        for (size_t row = 0; row < grid.rows(); ++row)
            for (size_t col = 0; col < grid.cols(); ++col)
                grid(row, col) = std::invoke(func, col * dx, row * dy);
    }

    template<typename Meta, typename Callable>
    constexpr void apply_equi2D(
        _basic_grid<Meta>& grid,
        const grid_value_t<Meta> x_min,
        const grid_value_t<Meta> x_max,
        Callable func)
    {
        return apply_equi2D(grid, x_min, x_max, x_min, x_max, func);
    }

}


#endif // !LLPS_GRID_HPP_INCLUDED