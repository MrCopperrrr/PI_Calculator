#pragma once

#ifdef CORE_EXPORTS
// This directive tells the compiler that the functions are being exported from this DLL.
// This is what causes the .lib file to be generated.
#define CORE_API __declspec(dllexport) 
#else
// This tells projects using the DLL that the functions are being imported.
#define CORE_API __declspec(dllimport)
#endif

#include <string>

namespace PiCalculator
{
    class Chudnovsky
    {
    public:
        // The CORE_API macro will be expanded to __declspec(dllexport) here
        static CORE_API std::string calculate(int digits, int num_threads);
    };
}