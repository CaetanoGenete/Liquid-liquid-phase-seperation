#ifndef LLPS_TIMER_HPP_INCLUDED
#define LLPS_TIMER_HPP_INCLUDED

#include <chrono>
#include <iostream>
#include <ctime>

namespace llps {

    struct timer
    {
    public:
        timer():
            _started_at(std::chrono::steady_clock::now()) 
        {
            auto start_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
            std::cout << "Computation started at: " << ctime(&start_time) << "\n";
        }

        ~timer()
        {
            auto finish_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
            std::cout << "Computation finited at: " << ctime(&finish_time) << "\n";

            std::chrono::duration<double> elapsed_time = std::chrono::steady_clock::now() - _started_at;
            std::chrono::hh_mm_ss formated_time{ elapsed_time };

            std::cout << "Time elapsed (hh:mm:ss): " << formated_time << "\n";

        }
    private:
        std::chrono::steady_clock::time_point _started_at;
    };

}

#endif // !LLPS_TIMER_HPP_INCLUDED