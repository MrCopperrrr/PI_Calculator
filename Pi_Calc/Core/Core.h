#pragma once

#ifdef CORE_EXPORTS
#define CORE_API __declspec(dlspecport)
#else
#define CORE_API __declspec(dllimport)
#endif

#include <string>

namespace PiCalculator
{
    class Chudnovsky
    {
    public:
        static std::string calculate(int digits, int num_threads);
    };
}