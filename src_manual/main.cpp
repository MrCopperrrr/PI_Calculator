#include "bigint.hpp"
#include "timer.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <omp.h>
#include <string>
#include <vector>


using namespace pi;

int64_t parse_digits(std::string arg) {
  if (arg.empty())
    return 1000;
  char suffix = std::tolower(arg.back());
  int64_t multiplier = 1;
  bool has_suffix = true;
  if (suffix == 'k')
    multiplier = 1000;
  else if (suffix == 'm')
    multiplier = 1000000;
  else if (suffix == 'b')
    multiplier = 1000000000;
  else
    has_suffix = false;
  std::string num_part = has_suffix ? arg.substr(0, arg.size() - 1) : arg;
  try {
    return std::stoll(num_part) * multiplier;
  } catch (...) {
    return 1000;
  }
}

int main(int argc, char *argv[]) {
  int64_t digits = 1000;
  if (argc > 1)
    digits = parse_digits(argv[1]);

  Timer total_timer;
  int64_t iterations = digits / 14 + 1;

  std::cout << "--- PI Calculator (Chudnovsky NTT) ---" << std::endl;
  std::cout << "Target: " << digits << " digits" << std::endl;

  BigInt P, Q, T;
  Timer step_timer;
  std::cout << "Step 1: Binary Splitting... ";
#pragma omp parallel
  {
#pragma omp single
    {
      BigInt::binary_split(0, iterations, P, Q, T);
    }
  }
  std::cout << "Done in " << step_timer.elapsed_str() << std::endl;

  // Pi = (426880 * sqrt(10005) * Q) / T
  // For now, we report the limbs.
  // Implementing Newton-Raphson Sqrt and Div...

  std::cout << "Step 2: Finalizing Result (Placeholder for NR)..." << std::endl;
  std::cout << "  - Q limbs: " << Q.data.size() << " (approx "
            << Q.data.size() * 19 << " digits)" << std::endl;
  std::cout << "  - T limbs: " << T.data.size() << std::endl;

  std::cout << "-----------------------------------" << std::endl;
  std::cout << "Total execution time: " << total_timer.elapsed_str()
            << std::endl;

  return 0;
}
