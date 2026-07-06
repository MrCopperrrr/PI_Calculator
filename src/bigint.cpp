#include "bigint.hpp"
#include "ntt.hpp"
#include <omp.h>
#include <vector>

namespace pi {

void parallel_reciprocal(mpz_t R, const mpz_t X, size_t k) {
  size_t x_bits = mpz_sizeinbase(X, 2);
  if (k < 2000 || x_bits < 2000) {
    mpz_t t;
    mpz_init_set_ui(t, 1);
    mpz_mul_2exp(t, t, k);
    mpz_tdiv_q(R, t, X);
    mpz_clear(t);
    return;
  }

  mpz_t X_small, R_small;
  mpz_inits(X_small, R_small, NULL);

  size_t half_k = k / 2 + 128;
  size_t shift = (x_bits > half_k) ? (x_bits - half_k) : 0;
  
  mpz_tdiv_q_2exp(X_small, X, shift);
  parallel_reciprocal(R_small, X_small, half_k);

  // Scale R_small back up
  // R \approx 2^k / X
  // R_small \approx 2^half_k / X_small = 2^half_k / (X >> shift) \approx 2^half_k / (X / 2^shift)
  // R_small \approx 2^(half_k + shift) / X
  // So R \approx R_small * 2^(k - half_k - shift)
  if (k > half_k + shift) {
    mpz_mul_2exp(R_small, R_small, k - half_k - shift);
  } else {
    mpz_tdiv_q_2exp(R_small, R_small, half_k + shift - k);
  }
  mpz_set(R, R_small);

  // Newton iteration: R = R + R(2^k - XR) / 2^k
  mpz_t E, T, p2k;
  mpz_inits(E, T, p2k, NULL);

  NTTMultiplier::multiply(E, X, R);
  mpz_set_ui(p2k, 1);
  mpz_mul_2exp(p2k, p2k, k);
  mpz_sub(E, p2k, E);

  NTTMultiplier::multiply(T, R, E);
  mpz_tdiv_q_2exp(T, T, k);
  mpz_add(R, R, T);

  mpz_clears(X_small, R_small, E, T, p2k, NULL);
}

void parallel_invsqrt(mpz_t R, const mpz_t X, size_t k) {
  size_t x_bits = mpz_sizeinbase(X, 2);
  if (k < 2000 || x_bits < 2000) {
    mpz_t t;
    mpz_init(t);
    mpz_sqrt(t, X);
    mpz_t r;
    mpz_init_set_ui(r, 1);
    mpz_mul_2exp(r, r, k);
    mpz_tdiv_q(R, r, t);
    mpz_clears(t, r, NULL);
    return;
  }

  mpz_t X_small, R_small;
  mpz_inits(X_small, R_small, NULL);

  size_t half_k = k / 2 + 128;
  size_t shift = (x_bits > 2 * half_k) ? (x_bits - 2 * half_k) : 0;
  if (shift % 2 != 0) shift++; 

  mpz_tdiv_q_2exp(X_small, X, shift);
  parallel_invsqrt(R_small, X_small, half_k);

  // R \approx 2^k / sqrt(X)
  // R_small \approx 2^half_k / sqrt(X_small) \approx 2^half_k / sqrt(X / 2^shift) = 2^(half_k + shift/2) / sqrt(X)
  size_t target_shift = half_k + shift / 2;
  if (k > target_shift) {
    mpz_mul_2exp(R_small, R_small, k - target_shift);
  } else {
    mpz_tdiv_q_2exp(R_small, R_small, target_shift - k);
  }
  mpz_set(R, R_small);

  // Newton iteration: R = R * (3*2^(2k) - X*R^2) / 2^(2k+1)
  mpz_t R2, XR2, p2k3, T;
  mpz_inits(R2, XR2, p2k3, T, NULL);

  NTTMultiplier::multiply(R2, R, R);
  NTTMultiplier::multiply(XR2, X, R2);
  
  mpz_set_ui(p2k3, 3);
  mpz_mul_2exp(p2k3, p2k3, 2 * k);
  mpz_sub(T, p2k3, XR2);
  
  NTTMultiplier::multiply(R2, R, T);
  mpz_tdiv_q_2exp(R, R2, 2 * k + 1);

  mpz_clears(X_small, R_small, R2, XR2, p2k3, T, NULL);
}

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
  size_t bits = mpz_sizeinbase(n, 2);
  if (bits < 1000000000) { // GMP assembly is extremely fast for roots < 300M digits
    mpz_sqrt(rop, n);
    return;
  }
  // ... rest of the function ...

  // To get a correct sqrt, the reciprocal needs to be scaled properly
  size_t k = bits + 64; 
  mpz_t inv;
  mpz_init(inv);
  
#pragma omp parallel
  {
#pragma omp single
    parallel_invsqrt(inv, n, k);
  }

  NTTMultiplier::multiply(rop, n, inv);
  mpz_tdiv_q_2exp(rop, rop, k);
  
  // Refine with GMP's native sqrt for final bits (very fast)
  mpz_t rem;
  mpz_init(rem);
  mpz_mul(rem, rop, rop);
  mpz_sub(rem, n, rem);
  if (mpz_sgn(rem) < 0) {
      mpz_sqrt(rop, n); // Fallback to safe version if Newton was off
  } else {
      // Newton is usually floor or floor-1
      mpz_t next_rop;
      mpz_init(next_rop);
      mpz_add_ui(next_rop, rop, 1);
      mpz_mul(next_rop, next_rop, next_rop);
      if (mpz_cmp(next_rop, n) <= 0) {
          mpz_sqrt(rop, n); // Just to be safe
      }
      mpz_clear(next_rop);
  }
  mpz_clear(rem);
  mpz_clear(inv);
}

void BigInt::parallel_div(mpz_t q, const mpz_t num, const mpz_t den) {
  size_t n_bits = mpz_sizeinbase(num, 2);
  size_t d_bits = mpz_sizeinbase(den, 2);
  
  if (n_bits < 1000000 || d_bits < 500000) {
    mpz_tdiv_q(q, num, den);
    return;
  }

  size_t k = n_bits + 64;
  mpz_t inv;
  mpz_init(inv);
  
#pragma omp parallel
  {
#pragma omp single
    parallel_reciprocal(inv, den, k);
  }

  NTTMultiplier::multiply(q, num, inv);
  mpz_tdiv_q_2exp(q, q, k);
  
  // Correction step: Newton gives a very close approximation
  mpz_t rem, prod;
  mpz_inits(rem, prod, NULL);
  mpz_mul(prod, q, den);
  mpz_sub(rem, num, prod);
  
  // If the remainder is still larger than divisor, or negative,
  // we use one final mpz_tdiv_q on the remainder for perfect accuracy.
  if (mpz_cmp(rem, den) >= 0 || mpz_sgn(rem) < 0) {
    mpz_t q_corr;
    mpz_init(q_corr);
    mpz_tdiv_q(q_corr, rem, den);
    mpz_add(q, q, q_corr);
    mpz_clear(q_corr);
  }
  
  mpz_clears(inv, rem, prod, NULL);
}

} // namespace pi
