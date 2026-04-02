#include "bigint.hpp"
#include "ntt.hpp"
#include <omp.h>
#include <vector>

namespace pi {

// Tree-based Parallel Power: 10^N = 10^(N/2) * 10^(N/2)
void recursive_pow(mpz_t rop, uint64_t base, uint64_t exp) {
  if (exp == 0) {
    mpz_set_ui(rop, 1);
    return;
  }
  if (exp == 1) {
    mpz_set_ui(rop, base);
    return;
  }

  mpz_t half;
  mpz_init(half);

  // For large exponents, split and calculate in parallel
  if (exp > 100000) {
    mpz_t res2;
    mpz_init(res2);

#pragma omp task shared(half)
    recursive_pow(half, base, exp / 2);

#pragma omp task shared(res2)
    recursive_pow(res2, base, exp - (exp / 2));

#pragma omp taskwait

    NTTMultiplier::multiply(rop, half, res2);
    mpz_clear(res2);
  } else {
    // Fallback to GMP's native power for smaller subproblems
    mpz_ui_pow_ui(rop, base, exp);
  }
  mpz_clear(half);
}

void BigInt::parallel_pow_ui(mpz_t rop, uint64_t base, uint64_t exp) {
#pragma omp parallel
  {
#pragma omp single
    recursive_pow(rop, base, exp);
  }
}

void BigInt::parallel_sqrt(mpz_t rop, const mpz_t n) {
  if (mpz_sgn(n) <= 0) {
    mpz_set_ui(rop, 0);
    return;
  }
  // GMP's sqrt is already very fast for 100M.
  // The main bottleneck is the multiplication after it.
  mpz_sqrt(rop, n);
}

} // namespace pi
