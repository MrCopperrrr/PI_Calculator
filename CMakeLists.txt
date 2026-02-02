cmake_minimum_required(VERSION 3.10)
project(PI_Calculator)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MSVC)
    add_compile_options(/O2 /arch:AVX2 /openmp)
else()
    add_compile_options(-O3 -mavx2 -march=native -ffast-math -pthread -fopenmp)
endif()

find_package(OpenMP REQUIRED)

set(SOURCES 
    src/main.cpp 
    src/bigint.cpp 
    src/chudnovsky.cpp 
    src/base_conv.cpp
)

include_directories(include)

find_library(GMP_LIB NAMES gmp libgmp)

add_executable(pi_calc ${SOURCES})
target_link_libraries(pi_calc PRIVATE ${GMP_LIB} OpenMP::OpenMP_CXX)
