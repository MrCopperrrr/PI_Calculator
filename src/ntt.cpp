#include "ntt.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <gmp.h>
#include <vector>

namespace pi {

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

void parallel_mul_recursive(mpz_t rop, const mpz_t op1, const mpz_t op2,
                            int depth) {
  size_t bits = std::max(mpz_sizeinbase(op1, 2), mpz_sizeinbase(op2, 2));

  // Only split if depth > 0 AND combined bit count is large enough to justify
  // overhead
  if (depth <= 0 || bits < 4000000) {
    mpz_mul(rop, op1, op2);
    return;
  }

  size_t split = bits / 2;
  mpz_t a_h, a_l, b_h, b_l, r1, r2, r3, r4;
  mpz_inits(a_h, a_l, b_h, b_l, r1, r2, r3, r4, NULL);

  mpz_tdiv_q_2exp(a_h, op1, split);
  mpz_tdiv_r_2exp(a_l, op1, split);
  mpz_tdiv_q_2exp(b_h, op2, split);
  mpz_tdiv_r_2exp(b_l, op2, split);

#pragma omp task shared(r1) firstprivate(depth)
  parallel_mul_recursive(r1, a_h, b_h, depth - 1);
#pragma omp task shared(r2) firstprivate(depth)
  parallel_mul_recursive(r2, a_h, b_l, depth - 1);
#pragma omp task shared(r3) firstprivate(depth)
  parallel_mul_recursive(r3, a_l, b_h, depth - 1);
#pragma omp task shared(r4) firstprivate(depth)
  parallel_mul_recursive(r4, a_l, b_l, depth - 1);

#pragma omp taskwait

  mpz_mul_2exp(r1, r1, 2 * split);
  mpz_add(r2, r2, r3);
  mpz_mul_2exp(r2, r2, split);
  mpz_add(rop, r1, r2);
  mpz_add(rop, rop, r4);

  mpz_clears(a_h, a_l, b_h, b_l, r1, r2, r3, r4, NULL);
}

void NTTMultiplier::multiply(mpz_t rop, const mpz_t op1, const mpz_t op2) {
  // Threshold to even enter parallel region: 4 Million bits
  if (std::max(mpz_sizeinbase(op1, 2), mpz_sizeinbase(op2, 2)) < 4000000) {
    mpz_mul(rop, op1, op2);
    return;
  }

#pragma omp parallel
  {
#pragma omp single
    parallel_mul_recursive(rop, op1, op2, 2);
  }
}

} // namespace pi
