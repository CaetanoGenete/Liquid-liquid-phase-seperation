    #ifndef LLPS_ARRAY_OF_RANGES_ALGEBRA_HPP_INCLUDED
    #define LLPS_ARRAY_OF_RANGES_ALGEBRA_HPP_INCLUDED

    #include <boost/numeric/odeint/algebra/range_algebra.hpp>

    struct array_of_ranges_algebra
    {
        using algebra = boost::numeric::odeint::range_algebra;        
        
        template<template<typename, size_t> class Array, typename T, size_t dim, class Op>
    static void for_each1(Array<T, dim>& s1, , Op op)
    {
        for(size_t i = 0; i < dim; ++i)
            algebra::for_each1(s1[i], op);
    }

    template<template<typename, size_t> class Array, typename T, size_t dim, class Op>
    static void for_each2(Array<T, dim>& s1, const Array<T, dim>& s2, Op op)
    {
        for(size_t i = 0; i < dim; ++i)
            algebra::for_each2(s1[i], s2[i], op);
    }

    /* different const signature - required for the scale_sum_swap2 operation */
    template<template<typename, size_t> class Array, typename T, size_t dim, class Op>
    static void for_each3(Array<T, dim>& s1, Array<T, dim>& s2, const Array<T, dim>& s3, Op op)
    {
        for(size_t i = 0; i < dim; ++i)
            algebra::for_each3(s1[i], s2[i], s3[i], op);
    }

    template<template<typename, size_t> class Array, typename T, size_t dim, class Op>
    static void for_each4(Array<T, dim>& s1, const Array<T, dim>& s2, const Array<T, dim>& s3, const Array<T, dim>& s4, Op op)
    {
        for(size_t i = 0; i < dim; ++i)
            algebra::for_each4(s1[i], s2[i], s3[i], s4[i], op);
    }

    template<template<typename, size_t> class Array, typename T, size_t dim, class Op>
    static void for_each5(Array<T, dim>& s1, const Array<T, dim>& s2, const Array<T, dim>& s3, const Array<T, dim>& s4, const Array<T, dim>& s5, Op op)
    {
        for(size_t i = 0; i < dim; ++i)
            algebra::for_each5(s1[i], s2[i], s3[i], s4[i], s5[i], op);
    }

    template<template<typename, size_t> class Array, typename T, size_t dim, class Op>
    static void for_each6(Array<T, dim>& s1, const Array<T, dim>& s2, const Array<T, dim>& s3, const Array<T, dim>& s4, const Array<T, dim>& s5, const Array<T, dim>& s6, Op op)
    {
        for(size_t i = 0; i < dim; ++i)
            algebra::for_each6(s1[i], s2[i], s3[i], s4[i], s5[i], s6[i], op);
    }

    template<template<typename, size_t> class Array, typename T, size_t dim, class Op>
    static void for_each7(Array<T, dim>& s1, const Array<T, dim>& s2, const Array<T, dim>& s3, const Array<T, dim>& s4, const Array<T, dim>& s5, const Array<T, dim>& s6, const Array<T, dim>& s7, Op op)
    {
        for(size_t i = 0; i < dim; ++i)
            algebra::for_each7(s1[i], s2[i], s3[i], s4[i], s5[i], s6[i], s7[i], op);
    }

    template<template<typename, size_t> class Array, typename T, size_t dim, class Op>
    static void for_each8(Array<T, dim>& s1, const Array<T, dim>& s2, const Array<T, dim>& s3, const Array<T, dim>& s4, const Array<T, dim>& s5, const Array<T, dim>& s6, const Array<T, dim>& s7, const Array<T, dim>& s8, Op op)
    {
        for(size_t i = 0; i < dim; ++i)
            algebra::for_each8(s1[i], s2[i], s3[i], s4[i], s5[i], s6[i], s7[i], s8[i], op);
    }

    template<template<typename, size_t> class Array, typename T, size_t dim, class Op>
    static void for_each9(Array<T, dim>& s1, const Array<T, dim>& s2, const Array<T, dim>& s3, const Array<T, dim>& s4, const Array<T, dim>& s5, const Array<T, dim>& s6, const Array<T, dim>& s7, const Array<T, dim>& s8, const Array<T, dim>& s9, Op op)
    {
        for(size_t i = 0; i < dim; ++i)
            algebra::for_each9(s1[i], s2[i], s3[i], s4[i], s5[i], s6[i], s7[i], s8[i], s9[i], op);
    }

    template<template<typename, size_t> class Array, typename T, size_t dim, class Op>
    static void for_each10(Array<T, dim>& s1, const Array<T, dim>& s2, const Array<T, dim>& s3, const Array<T, dim>& s4, const Array<T, dim>& s5, const Array<T, dim>& s6, const Array<T, dim>& s7, const Array<T, dim>& s8, const Array<T, dim>& s9, const Array<T, dim>& s10, Op op)
    {
        for(size_t i = 0; i < dim; ++i)
            algebra::for_each10(s1[i], s2[i], s3[i], s4[i], s5[i], s6[i], s7[i], s8[i], s9[i], s10[i], op);
    }

    template<template<typename, size_t> class Array, typename T, size_t dim, class Op>
    static void for_each11(Array<T, dim>& s1, const Array<T, dim>& s2, const Array<T, dim>& s3, const Array<T, dim>& s4, const Array<T, dim>& s5, const Array<T, dim>& s6, const Array<T, dim>& s7, const Array<T, dim>& s8, const Array<T, dim>& s9, const Array<T, dim>& s10, const Array<T, dim>& s11, Op op)
    {
        for(size_t i = 0; i < dim; ++i)
            algebra::for_each11(s1[i], s2[i], s3[i], s4[i], s5[i], s6[i], s7[i], s8[i], s9[i], s10[i], s11[i], op);
    }

    template<template<typename, size_t> class Array, typename T, size_t dim, class Op>
    static void for_each12(Array<T, dim>& s1, const Array<T, dim>& s2, const Array<T, dim>& s3, const Array<T, dim>& s4, const Array<T, dim>& s5, const Array<T, dim>& s6, const Array<T, dim>& s7, const Array<T, dim>& s8, const Array<T, dim>& s9, const Array<T, dim>& s10, const Array<T, dim>& s11, const Array<T, dim>& s12, Op op)
    {
        for(size_t i = 0; i < dim; ++i)
            algebra::for_each12(s1[i], s2[i], s3[i], s4[i], s5[i], s6[i], s7[i], s8[i], s9[i], s10[i], s11[i], s12[i], op);
    }

    template<template<typename, size_t> class Array, typename T, size_t dim, class Op>
    static void for_each13(Array<T, dim>& s1, const Array<T, dim>& s2, const Array<T, dim>& s3, const Array<T, dim>& s4, const Array<T, dim>& s5, const Array<T, dim>& s6, const Array<T, dim>& s7, const Array<T, dim>& s8, const Array<T, dim>& s9, const Array<T, dim>& s10, const Array<T, dim>& s11, const Array<T, dim>& s12, const Array<T, dim>& s13, Op op)
    {
        for(size_t i = 0; i < dim; ++i)
            algebra::for_each13(s1[i], s2[i], s3[i], s4[i], s5[i], s6[i], s7[i], s8[i], s9[i], s10[i], s11[i], s12[i], s13[i], op);
    }

    template<template<typename, size_t> class Array, typename T, size_t dim, class Op>
    static void for_each14(Array<T, dim>& s1, const Array<T, dim>& s2, const Array<T, dim>& s3, const Array<T, dim>& s4, const Array<T, dim>& s5, const Array<T, dim>& s6, const Array<T, dim>& s7, const Array<T, dim>& s8, const Array<T, dim>& s9, const Array<T, dim>& s10, const Array<T, dim>& s11, const Array<T, dim>& s12, const Array<T, dim>& s13, const Array<T, dim>& s14, Op op)
    {
        for(size_t i = 0; i < dim; ++i)
            algebra::for_each14(s1[i], s2[i], s3[i], s4[i], s5[i], s6[i], s7[i], s8[i], s9[i], s10[i], s11[i], s12[i], s13[i], s14[i], op);
    }

    template<template<typename, size_t> class Array, typename T, size_t dim, class Op>
    static void for_each15(Array<T, dim>& s1, const Array<T, dim>& s2, const Array<T, dim>& s3, const Array<T, dim>& s4, const Array<T, dim>& s5, const Array<T, dim>& s6, const Array<T, dim>& s7, const Array<T, dim>& s8, const Array<T, dim>& s9, const Array<T, dim>& s10, const Array<T, dim>& s11, const Array<T, dim>& s12, const Array<T, dim>& s13, const Array<T, dim>& s14, const Array<T, dim>& s15, Op op)
    {
        for(size_t i = 0; i < dim; ++i)
            algebra::for_each15(s1[i], s2[i], s3[i], s4[i], s5[i], s6[i], s7[i], s8[i], s9[i], s10[i], s11[i], s12[i], s13[i], s14[i], s15[i], op);
    }

};

#endif // !LLPS_ARRAY_OF_RANGES_ALGEBRA_HPP_INCLUDED