#include "ntt.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <gmp.h>
#include <vector>

namespace pi {
bool NTTMultiplier::use_hybrid = false;

uint64_t NTTMultiplier::power(uint64_t base, uint64_t exp, uint64_t mod) {
  uint64_t res = 1;
  base %= mod;
  while (exp > 0) {
    if (exp % 2 == 1)
      res = (__uint128_t)res * base % mod;
    base = (__uint128_t)base * base % mod;
    exp /= 2;
  }
  return res;
}

uint64_t NTTMultiplier::modInverse(uint64_t n, uint64_t mod) {
  return power(n, mod - 2, mod);
}

void NTTMultiplier::ntt(std::vector<uint64_t> &a, bool invert, uint64_t mod) {
  int n = a.size();
  for (int i = 1, j = 0; i < n; i++) {
    int bit = n >> 1;
    for (; j & bit; bit >>= 1)
      j ^= bit;
    j ^= bit;
    if (i < j)
      std::swap(a[i], a[j]);
  }

  for (int len = 2; len <= n; len <<= 1) {
    uint64_t wlen = power(G, (mod - 1) / len, mod);
    if (invert)
      wlen = modInverse(wlen, mod);

#pragma omp parallel for if (n > 65536) schedule(static)
    for (int i = 0; i < n; i += len) {
      uint64_t w = 1;
      for (int j = 0; j < len / 2; j++) {
        uint64_t u = a[i + j];
        uint64_t v = (__uint128_t)a[i + j + len / 2] * w % mod;
        a[i + j] = (u + v) % mod;
        if (a[i + j] >= mod)
          a[i + j] -= mod;
        a[i + j + len / 2] = (u - v + mod) % mod;
        w = (__uint128_t)w * wlen % mod;
      }
    }
  }

  if (invert) {
    uint64_t n_inv = modInverse(n, mod);
    for (uint64_t &x : a)
      x = (__uint128_t)x * n_inv % mod;
  }
}

void parallel_mul_karatsuba(mpz_t rop, const mpz_t op1, const mpz_t op2,
                            int depth) {
  size_t bits1 = mpz_sizeinbase(op1, 2);
  size_t bits2 = mpz_sizeinbase(op2, 2);
  size_t max_bits = std::max(bits1, bits2);

  if (depth <= 0 || max_bits < 4000000) {
    mpz_mul(rop, op1, op2);
    return;
  }

  size_t split = max_bits / 2;

  mpz_t a_h, a_l, b_h, b_l;
  mpz_inits(a_h, a_l, b_h, b_l, NULL);

  mpz_tdiv_q_2exp(a_h, op1, split);
  mpz_tdiv_r_2exp(a_l, op1, split);
  mpz_tdiv_q_2exp(b_h, op2, split);
  mpz_tdiv_r_2exp(b_l, op2, split);

  mpz_t z2, z0, z1, sum_a, sum_b;
  mpz_inits(z2, z0, z1, sum_a, sum_b, NULL);

  mpz_add(sum_a, a_h, a_l);
  mpz_add(sum_b, b_h, b_l);

#pragma omp task shared(z2)
  parallel_mul_karatsuba(z2, a_h, b_h, depth - 1);

#pragma omp task shared(z0)
  parallel_mul_karatsuba(z0, a_l, b_l, depth - 1);

#pragma omp task shared(z1)
  parallel_mul_karatsuba(z1, sum_a, sum_b, depth - 1);

#pragma omp taskwait

  mpz_sub(z1, z1, z2);
  mpz_sub(z1, z1, z0);

  mpz_mul_2exp(z2, z2, 2 * split);
  mpz_mul_2exp(z1, z1, split);
  mpz_add(rop, z2, z1);
  mpz_add(rop, rop, z0);

  mpz_clears(a_h, a_l, b_h, b_l, z2, z0, z1, sum_a, sum_b, NULL);
}

void NTTMultiplier::multiply(mpz_t rop, const mpz_t op1, const mpz_t op2) {
  size_t max_bits = std::max(mpz_sizeinbase(op1, 2), mpz_sizeinbase(op2, 2));

  if (mpz_sgn(op1) < 0 || mpz_sgn(op2) < 0 || max_bits < 500000) {
    mpz_mul(rop, op1, op2);
    return;
  }

  // Strategy for Step 1: Binary Splitting
  // For extremely large numbers (> 20M bits), GMP's native FFT (O(n log n)) 
  // is faster than our Parallel Karatsuba (O(n^1.58)) even on 1 core.
  if (!use_hybrid && max_bits > 20000000) {
    mpz_mul(rop, op1, op2);
    return;
  }

  // Strategy for Step 2: Evaluation (Newton-Raphson)
  // We MUST use parallel Karatsuba even for huge numbers to bypass GMP's 
  // single-threaded bottleneck during divisions/sqrt.
  int depth = 4;
  if (max_bits > 50000000) {
    depth = 3; // Limit depth to 3 levels (27 tasks) for huge numbers
  }

#pragma omp parallel
    {
#pragma omp single
    parallel_mul_recursive(rop, op1, op2, depth);
  }
}

} // namespace pi
