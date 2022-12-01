#ifndef LLPS_UTILITIES_IO_HPP_INCLUDED
#define LLPS_UTILITIES_IO_HPP_INCLUDED

#include <iterator> //For access to std::iter_value_t and iterator concepts
#include <fstream>  //For access to std::ofstream
#include <string>   //For access to std::string
#include <cassert>  //For access to assert macro
#include <cstddef>  //For access to fixed size types

namespace llps::utilities {

    template<typename Type>
    inline void serialise_to_binary(std::ofstream& stream, const Type& value)
    {
        stream.write(reinterpret_cast<const char*>(&value), sizeof(Type));
    }

    inline void serialise_string(std::ofstream& stream, std::string string)
    {
        serialise_to_binary<uint64_t>(stream, string.size());
        stream << string;
    }

    struct plot_header
    {
        std::string title   = "";
        std::string x_label = "";
        std::string y_label = "";
        std::string x_scale = "linear";
        std::string y_scale = "linear";
    };

    inline void serialise_plot_header(std::ofstream& stream, size_t count, plot_header meta = {})
    {
        serialise_string(stream, meta.title);
        serialise_string(stream, meta.x_label);
        serialise_string(stream, meta.y_label);
        serialise_string(stream, meta.x_scale);
        serialise_string(stream, meta.y_scale);

        serialise_to_binary<uint64_t>(stream, count);
    }

    struct line_header
    {
        std::string label     = "";
        std::string colour    = "#000000";
        std::string linestyle = "solid";
    };

    template<std::floating_point XType, std::floating_point YType>
    inline void serialise_line_header(std::ofstream& stream, size_t samples_count, line_header meta = {})
    {
        serialise_string(stream, meta.label);
        serialise_string(stream, meta.colour);
        serialise_string(stream, meta.linestyle);
        
        serialise_to_binary<uint8_t>(stream, sizeof(XType));
        serialise_to_binary<uint8_t>(stream, sizeof(YType));

        serialise_to_binary<uint64_t>(stream, samples_count);
    }

    template<std::floating_point ValueType, std::floating_point SpaceType>
    struct video_header
    {
        std::pair<ValueType, ValueType> min_max = {0., 0.};
        SpaceType dx = 1.;
        SpaceType dy = 1.;

        uint32_t interval = 20;
    };

    template<std::floating_point ValueType, std::floating_point SpaceType>
    inline void serialise_video_header(
        std::ofstream& stream,
        size_t width, size_t height,
        size_t frames,
        video_header<ValueType, SpaceType> meta = {})
    {
        serialise_to_binary<uint8_t>(stream, sizeof(ValueType));
        //serialise_to_binary<uint8_t>(stream, sizeof(TimeType));
        serialise_to_binary<uint8_t>(stream, sizeof(SpaceType));

        serialise_to_binary(stream, meta.min_max.first);
        serialise_to_binary(stream, meta.min_max.second);
        serialise_to_binary(stream, meta.dx);
        serialise_to_binary(stream, meta.dy);

        serialise_to_binary<uint64_t>(stream, width);
        serialise_to_binary<uint64_t>(stream, height);

        serialise_to_binary(stream, meta.interval);
        
        serialise_to_binary<uint64_t>(stream, frames);

    }
}

#endif // !LLPS_UTILITIES_IO_HPP_INCLUDED