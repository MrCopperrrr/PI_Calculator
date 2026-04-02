# Pi-Calc v2.5: High-Performance Parallel Computing Engine for Pi Calculation

## 1. Project Overview
Pi-Calc is a high-performance computational tool designed for the high-precision calculation of the mathematical constant Pi. By leveraging the Chudnovsky algorithm, the project implements advanced parallel computing techniques to scale from millions to billions of decimal digits. The engine is optimized for multi-core processors and efficient memory management, making it suitable for benchmarking hardware and measuring computational efficiency.

## 2. Core Architecture and Algorithms

### 2.1. Mathematical Foundation
The calculation is based on the Chudnovsky formula (1988), which provides approximately 14.18 digits of Pi per term. This series converges rapidly and is the standard for modern record-breaking Pi computations.

### 2.2. Binary Splitting Method
To handle the summation of the series efficiently, the project implements the Binary Splitting method. This approach transforms the sum into a product of large integers, reducing the overall computational complexity. 
- **Parallelization**: The binary splitting process is parallelized using OpenMP tasking, allowing recursive sub-tasks to be distributed across all available CPU cores.

### 2.3. Hybrid Multiplication Engine
The engine utilizes a custom hybrid multiplication strategy to bridge the gap between standard library performance and parallel requirements:
- **GMP Integration**: Leverages the GNU Multiple Precision Arithmetic Library (GMP) for low-level high-precision integer arithmetic.
- **Parallel Recursive Splitting**: For extremely large operands (typically exceeding 4 million bits), the system employs a parallel recursive strategy to distribute the workload, overcoming the single-threaded limitations of standard library multiplication.
- **NTT Foundation**: Includes a Number Theoretic Transform (NTT) implementation as a framework for fast-convolution extensions.

### 2.4. Parallel Base Conversion
Binary-to-decimal conversion is often a bottleneck in high-precision calculations. Pi-Calc utilizes a parallel recursive division strategy based on powers of 10 to ensure that output generation scales linearly with data size.

## 3. Technical Specifications

- **Language**: C++17
- **Parallel Computing**: OpenMP 4.5+ (Task-based parallelism)
- **Arbitrary Precision**: GNU Multiple Precision Arithmetic Library (GMP)
- **Hardware Acceleration**: AVX2 instruction set optimizations
- **Memory Management**: Intelligent limb management to handle multi-gigabyte integer objects.

## 4. Performance Metrics and Benchmarking

The application provides detailed analytical reporting including:
- **Event Logging**: Real-time timestamps for each computational stage.
- **CPU Utilization**: Measurement of aggregate core usage and thread efficiency.
- **Wall Time vs. CPU Time**: Distinction between raw computation and I/O-bound operations (conversion and writing).

### Platform Considerations
- **Linux/WSL2**: Recommended for large-scale calculations (1B+ digits) due to 64-bit limb management.
- **Windows (MinGW-w64)**: Optimized for native execution with support for calculations up to 500 million digits.

## 5. Installation and Build Instructions

### Prerequisites
- CMake 3.10 or higher
- GCC 9+ or Clang 10+
- GMP Library (libgmp-dev)
- OpenMP-capable compiler

### Build Process (Linux / WSL2)
```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

### Build Process (Windows - MinGW)
```powershell
mkdir build
cd build
cmake .. -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release
cmake --build .
```

## 6. Usage

```bash
# Calculate 100 million digits
./pi_calc 100M

# Calculate 1 billion digits
./pi_calc 1B
```
The result is exported to `pi.txt` in the execution directory.

## 7. License
This project is licensed under the MIT License.
