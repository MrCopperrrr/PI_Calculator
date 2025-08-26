#pragma once

#ifdef GPU_MODULE_EXPORTS
#define GPU_API __declspec(dllexport)
#else
#define GPU_API __declspec(dllimport)
#endif

/**
 * @brief Initializes the CUDA context and verifies a compatible GPU is available.
 * @return True if initialization is successful, false otherwise.
 */
extern "C" GPU_API bool initialize_gpu_context();

/**
 * @brief Multiplies two large numbers on the GPU using FFT.
 * @param num1 Pointer to the first number as an array of unsigned ints.
 * @param size1 The number of elements in the first array.
 * @param num2 Pointer to the second number.
 * @param size2 The number of elements in the second array.
 * @param result Pointer to the array where the result will be stored.
 * @param result_size Input: The allocated size of the result array. Output: The actual size of the result.
 */
extern "C" GPU_API void multiply_on_gpu(
    const unsigned int* num1, int size1,
    const unsigned int* num2, int size2,
    unsigned int* result, int& result_size
);