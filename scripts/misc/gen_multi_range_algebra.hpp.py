import sys

if __name__ == "__main__":
    if len(sys.argv) > 2:
        raise ValueError("Expected atmost one argument!")

    out_path = "array_of_ranges_algebra.hpp"
    if len(sys.argv) == 2:
        out_path = sys.argv[1] + "/" + out_path

    indent = 4 * " "
    result = "#ifndef LLPS_ARRAY_OF_RANGES_ALGEBRA_HPP_INCLUDED\n"\
             "#define LLPS_ARRAY_OF_RANGES_ALGEBRA_HPP_INCLUDED\n"\
             "\n"\
             "#include <boost/numeric/odeint/algebra/range_algebra.hpp>\n"\
             "\n"\
             "struct array_of_ranges_algebra\n"\
             "{\n"+\
             indent + "using algebra = boost::numeric::odeint::range_algebra;\n\n"

    for_each_counts = 15
    for i in range(1, for_each_counts+1):
        if(i == 3):
            result += indent + "/* different const signature - required for the scale_sum_swap2 operation */\n"
        result += indent +  "template<template<typename, size_t> class Array, typename T, size_t dim, class Op>\n"
        result += indent + f"static void for_each{i}(Array<T, dim>& s1, "

        start=2
        if i == 3:
            result += "Array<T, dim>& s2, "
            start=3

        result += ", ".join([f"const Array<T, dim>& s{k}" for k in range(start, i+1)]) + ", Op op)\n"
        result += indent + "{\n"+\
                  2*indent +  "for(size_t i = 0; i < dim; ++i)\n"+\
                  2*indent + f"    algebra::for_each{i}(" + ", ".join([f"s{k+1}[i]" for k in range(i)]) + ", op);\n"+\
                  indent + "}\n\n"

    result += "};\n\n"
    result += "#endif // !LLPS_ARRAY_OF_RANGES_ALGEBRA_HPP_INCLUDED"

    with open(out_path, "w") as file:
        file.write(result)

