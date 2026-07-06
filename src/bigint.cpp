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

void parallel_reciprocal(mpz_t inv, const mpz_t B, size_t prec_bits) {
  size_t m = mpz_sizeinbase(B, 2);
  if (prec_bits <= 1500000) {
    mpz_t B_top;
    mpz_init(B_top);
    if (m >= prec_bits) {
      mpz_tdiv_q_2exp(B_top, B, m - prec_bits);
    } else {
      mpz_mul_2exp(B_top, B, prec_bits - m);
    }
    mpz_t num;
    mpz_init(num);
    mpz_set_ui(num, 1);
    mpz_mul_2exp(num, num, 2 * prec_bits);
    mpz_tdiv_q(inv, num, B_top);
    mpz_clears(B_top, num, NULL);
    return;
  }

  size_t half = (prec_bits + 1) / 2;
  mpz_t inv_half;
  mpz_init(inv_half);
  parallel_reciprocal(inv_half, B, half);

  size_t guard = 64;
  size_t need_b_bits = prec_bits + guard;
  mpz_t B_top;
  mpz_init(B_top);
  size_t shift_b = 0;
  if (m > need_b_bits) {
    shift_b = m - need_b_bits;
    mpz_tdiv_q_2exp(B_top, B, shift_b);
  } else {
    mpz_set(B_top, B);
  }
  size_t m_top = mpz_sizeinbase(B_top, 2);

  mpz_t BX, two_pow, term, prod;
  mpz_inits(BX, two_pow, term, prod, NULL);

  NTTMultiplier::multiply(BX, B_top, inv_half);
  size_t target_pow = m_top + half;
  mpz_set_ui(two_pow, 1);
  mpz_mul_2exp(two_pow, two_pow, target_pow + 1);

  mpz_sub(term, two_pow, BX);
  NTTMultiplier::multiply(prod, inv_half, term);

  size_t right_shift = m_top + 2 * half - prec_bits;
  mpz_tdiv_q_2exp(inv, prod, right_shift);

  mpz_clears(inv_half, B_top, BX, two_pow, term, prod, NULL);
}

void BigInt::parallel_div(mpz_t q, const mpz_t num, const mpz_t den) {
  size_t n_bits = mpz_sizeinbase(num, 2);
  size_t d_bits = mpz_sizeinbase(den, 2);

  if (mpz_cmp(num, den) < 0) {
    mpz_set_ui(q, 0);
    return;
  }

  if (d_bits < 1500000 || (n_bits - d_bits) < 1500000) {
    mpz_tdiv_q(q, num, den);
    return;
  }

  size_t k = n_bits - d_bits + 5;
  mpz_t inv;
  mpz_init(inv);
  parallel_reciprocal(inv, den, k);

  mpz_t prod, rem;
  mpz_inits(prod, rem, NULL);

  NTTMultiplier::multiply(prod, num, inv);
  mpz_tdiv_q_2exp(q, prod, d_bits + k);

  NTTMultiplier::multiply(prod, q, den);
  mpz_sub(rem, num, prod);

  while (mpz_sgn(rem) < 0) {
    mpz_sub_ui(q, q, 1);
    mpz_add(rem, rem, den);
  }
  while (mpz_cmp(rem, den) >= 0) {
    mpz_add_ui(q, q, 1);
    mpz_sub(rem, rem, den);
  }

  mpz_clears(inv, prod, rem, NULL);
}

void BigInt::parallel_sqrt(mpz_t rop, const mpz_t n) {
  if (mpz_sgn(n) <= 0) {
    mpz_set_ui(rop, 0);
    return;
  }

  size_t bits = mpz_sizeinbase(n, 2);
  if (bits < 3000000) {
    mpz_sqrt(rop, n);
    return;
  }

  size_t target_bits = (bits + 1) / 2;
  size_t shift = bits - target_bits;
  if (shift % 2 != 0)
    shift++;

  mpz_t n_top, s0;
  mpz_inits(n_top, s0, NULL);
  mpz_tdiv_q_2exp(n_top, n, shift);

  parallel_sqrt(s0, n_top);

  mpz_mul_2exp(s0, s0, shift / 2);

  mpz_t q, sum;
  mpz_inits(q, sum, NULL);

  parallel_div(q, n, s0);
  mpz_add(sum, s0, q);
  mpz_tdiv_q_2exp(rop, sum, 1);

  mpz_t sq, next_sq;
  mpz_inits(sq, next_sq, NULL);
  NTTMultiplier::multiply(sq, rop, rop);

  while (mpz_cmp(sq, n) > 0) {
    mpz_sub_ui(rop, rop, 1);
    NTTMultiplier::multiply(sq, rop, rop);
  }

  mpz_add(next_sq, sq, rop);
  mpz_add(next_sq, next_sq, rop);
  mpz_add_ui(next_sq, next_sq, 1);
  while (mpz_cmp(next_sq, n) <= 0) {
    mpz_add_ui(rop, rop, 1);
    mpz_set(sq, next_sq);
    mpz_add(next_sq, sq, rop);
    mpz_add(next_sq, next_sq, rop);
    mpz_add_ui(next_sq, next_sq, 1);
  }

  mpz_clears(n_top, s0, q, sum, sq, next_sq, NULL);
}

} // namespace pi
