#include "gpu_multiplier.h"

#include <cufft.h>
#include <cuda_runtime.h>
#include <stdio.h> // For error printing with fprintf

// --- GPU INITIALIZATION ---

extern "C" GPU_API bool initialize_gpu_context()
{
    // A safe and common way to initialize the CUDA context is to call a lightweight function.
    cudaError_t err = cudaFree(0);

    if (err != cudaSuccess) {
        // If it fails, print a descriptive error to the console
        fprintf(stderr, "CUDA Error: Failed to initialize context: %s\n", cudaGetErrorString(err));
        return false;
    }

    // You can also add a check for the device count
    int deviceCount = 0;
    err = cudaGetDeviceCount(&deviceCount);
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA Error during device count: %s\n", cudaGetErrorString(err));
        return false;
    }

    if (deviceCount == 0) {
        fprintf(stderr, "CUDA Error: No CUDA-capable devices were found.\n");
        return false;
    }

    return true;
}

// --- GPU MULTIPLICATION KERNELS ---

// CUDA kernel to convert our uint array into a complex format for FFT and pad it with zeros.
__global__ void setup_transform_kernel(const unsigned int* input, cufftComplex* output, int input_size, int transform_size) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < transform_size) {
        if (idx < input_size) {
            output[idx].x = (float)input[idx];
            output[idx].y = 0.0f;
        }
        else {
            // Pad with zeros
            output[idx].x = 0.0f;
            output[idx].y = 0.0f;
        }
    }
}

// CUDA kernel to perform point-wise multiplication of the transformed numbers
__global__ void complex_pointwise_mult_kernel(cufftComplex* a, const cufftComplex* b, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        float ax = a[idx].x;
        float ay = a[idx].y;
        float bx = b[idx].x;
        float by = b[idx].y;
        a[idx].x = ax * bx - ay * by;
        a[idx].y = ax * by + ay * bx;
    }
}

// CUDA kernel to normalize the result of the inverse FFT and handle the "carries"
__global__ void normalize_and_carry_kernel(const cufftComplex* input, unsigned int* output, int transform_size) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < transform_size) {
        // Normalize by dividing by the transform size
        unsigned long long val = roundf(input[idx].x / transform_size);

        // Atomically add the value to the output array. This handles carries within a number base of 2^64.
        if (val > 0) {
            int current_idx = idx;
            unsigned long long carry = 0;

            // Add the current value to the first limb
            unsigned long long current_val = atomicAdd((unsigned long long*) & output[current_idx], val);
            carry = (current_val + val) / 0x100000000ULL; // Carry for base 2^32

            current_idx++;
            // Propagate the carry
            while (carry > 0 && current_idx < transform_size) {
                current_val = atomicAdd((unsigned long long*) & output[current_idx], carry);
                carry = (current_val + carry) / 0x100000000ULL;
                current_idx++;
            }
        }
    }
}


// --- MAIN GPU MULTIPLICATION FUNCTION ---

extern "C" GPU_API void multiply_on_gpu(
    const unsigned int* num1_host, int size1,
    const unsigned int* num2_host, int size2,
    unsigned int* result_host, int& result_size)
{
    if (size1 == 0 || size2 == 0) {
        result_size = 0;
        return;
    }

    int transform_size = 1;
    while (transform_size < size1 + size2) {
        transform_size *= 2;
    }

    cufftComplex* num1_device, * num2_device;
    cudaMalloc((void**)&num1_device, sizeof(cufftComplex) * transform_size);
    cudaMalloc((void**)&num2_device, sizeof(cufftComplex) * transform_size);

    unsigned int* result_device_uint;
    cudaMalloc((void**)&result_device_uint, sizeof(unsigned int) * (transform_size + 1)); // Extra limb for carry
    cudaMemset(result_device_uint, 0, sizeof(unsigned int) * (transform_size + 1));

    int threads_per_block = 256;
    int blocks_per_grid = (transform_size + threads_per_block - 1) / threads_per_block;

    setup_transform_kernel << <blocks_per_grid, threads_per_block >> > (num1_host, num1_device, size1, transform_size);
    setup_transform_kernel << <blocks_per_grid, threads_per_block >> > (num2_host, num2_device, size2, transform_size);

    cufftHandle plan;
    cufftPlan1d(&plan, transform_size, CUFFT_C2C, 1);

    cufftExecC2C(plan, num1_device, num1_device, CUFFT_FORWARD);
    cufftExecC2C(plan, num2_device, num2_device, CUFFT_FORWARD);

    complex_pointwise_mult_kernel << <blocks_per_grid, threads_per_block >> > (num1_device, num2_device, transform_size);

    cufftExecC2C(plan, num1_device, num1_device, CUFFT_INVERSE);

    normalize_and_carry_kernel << <blocks_per_grid, threads_per_block >> > (num1_device, result_device_uint, transform_size);

    cudaMemcpy(result_host, result_device_uint, sizeof(unsigned int) * (transform_size + 1), cudaMemcpyDeviceToHost);

    int final_size = transform_size + 1;
    while (final_size > 1 && result_host[final_size - 1] == 0) {
        final_size--;
    }
    result_size = final_size;

    cudaFree(num1_device);
    cudaFree(num2_device);
    cudaFree(result_device_uint);
    cufftDestroy(plan);
}