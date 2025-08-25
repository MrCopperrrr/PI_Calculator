#include "pqt.h"

static const long long TASK_THRESHOLD = 1024; 

PQT binary_split_parallel(long long a, long long b, const mpz_class &C3_over_24) {
    PQT res;
    if (b - a == 1) {
        long long k = a;
        if (k == 0) {
            res.P = mpz_class(1);
            res.Q = mpz_class(1);
            res.T = mpz_class(13591409); // A0
            return res;
        }
        // p(k) = (6k-5)(2k-1)(6k-1)
        mpz_class pk = mpz_class(static_cast<long>(6*k - 5));
        pk *= mpz_class(static_cast<long>(2*k - 1));
        pk *= mpz_class(static_cast<long>(6*k - 1));

        // q(k) = k^3 * (C^3 / 24)
        mpz_class kk = mpz_class(static_cast<long>(k));
        mpz_class qk = kk * kk * kk;
        qk *= C3_over_24;

        // a(k) = 13591409 + 545140134 * k
        mpz_class ak = mpz_class(static_cast<long>(13591409)) +
                       mpz_class(static_cast<long>(545140134)) * mpz_class(static_cast<long>(k));

        mpz_class tk = pk * ak;
        if (k & 1) tk = -tk;

        res.P = pk;
        res.Q = qk;
        res.T = tk;
        return res;
    }

    long long m = (a + b) >> 1;
    PQT L, R;

    bool spawn = (b - a) >= TASK_THRESHOLD;
    if (spawn) {
        #pragma omp task shared(L, C3_over_24)
        { L = binary_split_parallel(a, m, C3_over_24); }
        #pragma omp task shared(R, C3_over_24)
        { R = binary_split_parallel(m, b, C3_over_24); }
        #pragma omp taskwait
    } else {
        L = binary_split_parallel(a, m, C3_over_24);
        R = binary_split_parallel(m, b, C3_over_24);
    }

    combine(res, L, R);
    return res;
}
