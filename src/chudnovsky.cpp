#include "bigint.hpp"
#include "ntt.hpp"
#include <cstdint>
#include <gmp.h>
#include <iostream>
#include <omp.h>


namespace pi {

const uint64_t A = 13591409;
const uint64_t B = 545140134;

void calculate_base_single(int64_t k, mpz_t P, mpz_t Q, mpz_t T, mpz_t C3_24) {
  if (k == 0) {
    mpz_set_ui(P, 1);
    mpz_set_ui(Q, 1);
    mpz_set_ui(T, A);
  } else {
    mpz_set_ui(P, 6 * k - 5);
    mpz_mul_ui(P, P, 2 * k - 1);
    mpz_mul_ui(P, P, 6 * k - 1);
    mpz_neg(P, P);

    mpz_set_ui(Q, k);
    mpz_mul_ui(Q, Q, k);
    mpz_mul_ui(Q, Q, k);
    mpz_mul(Q, Q, C3_24);

    mpz_set_ui(T, B);
    mpz_mul_ui(T, T, k);
    mpz_add_ui(T, T, A);
    mpz_mul(T, T, P);
  }
}

static int64_t total_it = 0;
static int64_t completed_it = 0;

void BigInt::binary_split(int64_t a, int64_t b, BigInt &P, BigInt &Q,
                          BigInt &T) {
  static mpz_t C3_24;
  static bool init_done = false;
#pragma omp critical
  {
    if (!init_done) {
      mpz_init_set_str(C3_24, "10939058860032000", 10);
      total_it = b;
      completed_it = 0;
      init_done = true;
    }
  }

  if (b - a == 1) {
    calculate_base_single(a, P.value, Q.value, T.value, C3_24);
#pragma omp atomic
    completed_it++;
    if (completed_it % 1000000 == 0) {
#pragma omp critical
      {
        printf("\rStep 1 Progress: %lld / %lld terms", (long long)completed_it,
               (long long)total_it);
        fflush(stdout);
      }
    }
    return;
  }

  int64_t m = (a + b) / 2;

  // Optimized tasking: only spawn tasks for large enough chunks to save RAM
  if (b - a > 100000) {
    BigInt P1, Q1, T1, P2, Q2, T2;
#pragma omp task shared(P1, Q1, T1)
    binary_split(a, m, P1, Q1, T1);
#pragma omp task shared(P2, Q2, T2)
    binary_split(m, b, P2, Q2, T2);
#pragma omp taskwait

    mpz_t T_part2;
    mpz_init(T_part2);

    // High-level merge: use tasking instead of nested parallel regions
    // T = T1*Q2 + P1*T2, P = P1*P2, Q = Q1*Q2
#pragma omp task shared(T, T1, Q2)
    NTTMultiplier::multiply(T.value, T1.value, Q2.value);
#pragma omp task shared(T_part2, P1, T2)
    NTTMultiplier::multiply(T_part2, P1.value, T2.value);
#pragma omp task shared(P, P1, P2)
    NTTMultiplier::multiply(P.value, P1.value, P2.value);
#pragma omp task shared(Q, Q1, Q2)
    NTTMultiplier::multiply(Q.value, Q1.value, Q2.value);
#pragma omp taskwait

    // Early clear: these are no longer needed after the merge
    T1.clear();
    T2.clear();
    Q1.clear();
    P1.clear();
    P2.clear();
    Q2.clear();

    mpz_add(T.value, T.value, T_part2);
    mpz_clear(T_part2);
  } else {
    BigInt P1, Q1, T1, P2, Q2, T2;
    binary_split(a, m, P1, Q1, T1);
    binary_split(m, b, P2, Q2, T2);

    mpz_t T_part2;
    mpz_init(T_part2);
    
    // T = T1*Q2 + P1*T2
    mpz_mul(T.value, T1.value, Q2.value);
    mpz_mul(T_part2, P1.value, T2.value);
    mpz_add(T.value, T.value, T_part2);
    mpz_clear(T_part2);
    
    // Clear T1 and T2 immediately
    T1.clear();
    T2.clear();

    // P = P1*P2
    mpz_mul(P.value, P1.value, P2.value);
    P1.clear();
    P2.clear();

    // Q = Q1*Q2
    mpz_mul(Q.value, Q1.value, Q2.value);
    Q1.clear();
    Q2.clear();
  }
}

} // namespace pi
