#pragma once

#ifdef CORE_EXPORTS
#define CORE_API __declspec(dllexport)
#else
#define CORE_API __declspec(dllimport)
#endif

// The primary function exported by the Core.dll
extern "C" CORE_API void compute_pi(long long digits, int nthreads, const char* outfile);