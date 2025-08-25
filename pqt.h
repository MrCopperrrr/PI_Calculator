#include <bits/stdc++.h>
#include <gmp.h>
#include <gmpxx.h>
#include <mpfr.h>
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;
using namespace std::chrono;

struct PQT {
    mpz_class P;
    mpz_class Q;
    mpz_class T;
    PQT() : P(0), Q(0), T(0) {}
    PQT(const mpz_class &p, const mpz_class &q, const mpz_class &t) : P(p), Q(q), T(t) {}
};

// combine two PQT results: out = L ⊗ R
static inline void combine(PQT &out, const PQT &L, const PQT &R) {
    // P = P_L * P_R
    // Q = Q_L * Q_R
    // T = T_L * Q_R + P_L * T_R
    out.P = L.P * R.P;
    out.Q = L.Q * R.Q;
    out.T = L.T * R.Q + L.P * R.T;
}
