#pragma once
#include <gmp.h>
#include <gmpxx.h>

struct PQT {
    mpz_class P;
    mpz_class Q;
    mpz_class T;
    PQT() : P(0), Q(0), T(0) {}
    PQT(const mpz_class& p, const mpz_class& q, const mpz_class& t) : P(p), Q(q), T(t) {}
};

static inline void combine(PQT& out, const PQT& L, const PQT& R) {
    out.P = L.P * R.P;
    out.Q = L.Q * R.Q;
    out.T = L.T * R.Q + L.P * R.T;
}