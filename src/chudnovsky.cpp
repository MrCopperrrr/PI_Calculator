#include "bigint.hpp"
#include <cstdint>
#include <gmp.h>
#include <omp.h>


namespace pi {

const uint64_t A = 13591409;
const uint64_t B = 545140134;

void calculate_base(int64_t k, mpz_t P, mpz_t Q, mpz_t T, mpz_t C3_24) {
  if (k == 0) {
    mpz_set_ui(P, 1);
    mpz_set_ui(Q, 1);
    mpz_set_ui(T, A);
  } else {
    // P = -(6k-5)(2k-1)(6k-1)
    // Use mpz_set_ui and mpz_mul_ui to prevent 64-bit overflow
    mpz_set_ui(P, 6 * k - 5);
    mpz_mul_ui(P, P, 2 * k - 1);
    mpz_mul_ui(P, P, 6 * k - 1);
    mpz_neg(P, P);

    // Q = k^3 * (C^3 / 24)
    mpz_set_ui(Q, k);
    mpz_mul_ui(Q, Q, k);
    mpz_mul_ui(Q, Q, k);
    mpz_mul(Q, Q, C3_24);

    // T = (A + Bk) * P
    mpz_set_ui(T, B);
    mpz_mul_ui(T, T, k);
    mpz_add_ui(T, T, A);
    mpz_mul(T, T, P);
  }
}

void BigInt::binary_split(int64_t a, int64_t b, BigInt &P, BigInt &Q,
                          BigInt &T) {
  static mpz_t C3_24;
  static bool init = false;
#pragma omp critical
  {
    if (!init) {
      mpz_init_set_str(C3_24, "10939058860032000", 10);
      init = true;
    }
  }

  if (b - a == 1) {
    calculate_base(a, P.value, Q.value, T.value, C3_24);
    return;
  }

  int64_t m = (a + b) / 2;
  BigInt P1, Q1, T1, P2, Q2, T2;

  if (b - a > 8192) {
#pragma omp task shared(P1, Q1, T1)
    binary_split(a, m, P1, Q1, T1);
#pragma omp task shared(P2, Q2, T2)
    binary_split(m, b, P2, Q2, T2);
#pragma omp taskwait
  } else {
    binary_split(a, m, P1, Q1, T1);
    binary_split(m, b, P2, Q2, T2);
  }

  mpz_t T_part2;
  mpz_init(T_part2);

#pragma omp parallel sections if (b - a > 32768)
  {
#pragma omp section
    mpz_mul(T.value, T1.value, Q2.value);
#pragma omp section
    mpz_mul(T_part2, P1.value, T2.value);
#pragma omp section
    mpz_mul(P.value, P1.value, P2.value);
#pragma omp section
    mpz_mul(Q.value, Q1.value, Q2.value);
  }

  mpz_add(T.value, T.value, T_part2);
  mpz_clear(T_part2);
}

} // namespace pi
