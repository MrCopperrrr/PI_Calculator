#include "bigint.hpp"

namespace pi {

const uint64_t A = 13591409;
const uint64_t B = 545140134;
const uint64_t C = 640320;
const uint64_t C3_24 = 10939058860032000ULL;

void BigInt::binary_split(int64_t a, int64_t b, BigInt &P, BigInt &Q,
                          BigInt &T) {
  if (b - a == 1) {
    if (a == 0) {
      P = BigInt(1);
      Q = BigInt(1);
      T = BigInt(A);
    } else {
      uint64_t k = a;
      P = BigInt(6 * k - 5);
      P.mul_small(2 * k - 1);
      P.mul_small(6 * k - 1);
      P.negative = true; // P(k) = -(6k-5)(2k-1)(6k-1)

      Q = BigInt(k);
      Q.mul_small(k);
      Q.mul_small(k);
      Q.mul_small(C3_24);

      T = BigInt(A + B * k);
      T.mul(P);
    }
    return;
  }

  int64_t m = (a + b) / 2;
  BigInt P1, Q1, T1, P2, Q2, T2;

  if (b - a > 500) { // Threshold for parallel task
#pragma omp task shared(P1, Q1, T1)
    binary_split(a, m, P1, Q1, T1);
#pragma omp task shared(P2, Q2, T2)
    binary_split(m, b, P2, Q2, T2);
#pragma omp taskwait
  } else {
    binary_split(a, m, P1, Q1, T1);
    binary_split(m, b, P2, Q2, T2);
  }

  // T = T1 * Q2 + P1 * T2
  BigInt T1Q2 = T1;
  T1Q2.mul(Q2);

  BigInt P1T2 = P1;
  P1T2.mul(T2);

  T = T1Q2;
  T.add(P1T2);

  P = P1;
  P.mul(P2);
  Q = Q1;
  Q.mul(Q2);
}

} // namespace pi
