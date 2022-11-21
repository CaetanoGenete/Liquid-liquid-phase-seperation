#ifndef UTILITIES_IO_HPP_INCLUDED
#define UTILITIES_IO_HPP_INCLUDED

#include <iterator> //For access to std::iter_value_t and iterator concepts
#include <fstream>  //For access to std::ofstream
#include <string>   //For access to std::string
#include <cassert>  //For access to assert macro

namespace utilities {

    template<typename Type>
    inline void serialise_to_binary(std::ofstream& stream, const Type& value)
    {
        stream.write(reinterpret_cast<const char*>(&value), sizeof(Type));
    }

    inline void serialise_string(std::ofstream& stream, std::string string)
    {
        serialise_to_binary(stream, string.size());
        stream << string;
    }

    struct plot_meta
    {
        std::string title   = "";
        std::string x_label = "";
        std::string y_label = "";
    };

    inline void serialise_plot_header(std::ofstream& stream, size_t count, plot_meta meta = {})
    {
        serialise_string(stream, meta.title);
        serialise_string(stream, meta.x_label);
        serialise_string(stream, meta.y_label);

        serialise_to_binary(stream, count);
    }

    struct line_meta
    {
        std::string label     = "";
        std::string colour    = "#000000";
        std::string linestyle = "solid";
    };

    template<std::floating_point XType, std::floating_point YType>
    inline void serialise_line_header(std::ofstream& stream, size_t samples_count, line_meta meta = {})
    {
        serialise_string(stream, meta.label);
        serialise_string(stream, meta.colour);
        serialise_string(stream, meta.linestyle);
        
        serialise_to_binary<uint8_t>(stream, sizeof(XType));
        serialise_to_binary<uint8_t>(stream, sizeof(YType));

        serialise_to_binary(stream, samples_count);
    }
}

#endif // !UTILITIES_IO_HPP_INCLUDED