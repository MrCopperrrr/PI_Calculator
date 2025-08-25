#include "core_api.h"
#include "binary_split.h"
#include <mpfr.h>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector> // Needed for storing thread results
#include <algorithm> // For std::min

#ifdef _OPENMP
#include <omp.h>
#endif

// This is the SERIAL implementation of the binary split algorithm.
// All the OpenMP pragmas have been removed from this function.
PQT binary_split(long long a, long long b, const mpz_class& C3_over_24) {
    PQT res;
    if (b - a == 1) {
        long long k = a;
        if (k == 0) {
            res.P = 1;
            res.Q = 1;
            res.T = 13591409;
            return res;
        }
        // Cast k to long to resolve ambiguity with mpz_class operators
        mpz_class pk = (6 * static_cast<long>(k) - 5);
        pk *= (2 * static_cast<long>(k) - 1);
        pk *= (6 * static_cast<long>(k) - 1);

        mpz_class kk(static_cast<long>(k));
        mpz_class qk = kk * kk * kk;
        qk *= C3_over_24;

        mpz_class ak = 13591409 + mpz_class(545140134) * static_cast<long>(k);
        mpz_class tk = pk * ak;
        if (k & 1) tk = -tk;

        res.P = pk;
        res.Q = qk;
        res.T = tk;
        return res;
    }

    long long m = (a + b) / 2;
    PQT L = binary_split(a, m, C3_over_24);
    PQT R = binary_split(m, b, C3_over_24);

    combine(res, L, R);
    return res;
}

// Main exported function
void compute_pi(long long digits, int nthreads, const char* outfile) {
    if (digits < 1) digits = 1;
    if (nthreads < 1) nthreads = 1;

    long long N = (digits + 13) / 14;
    long long prec_bits = (long long)ceil(digits * 3.3219280948873626);
    long long extra_guard = 128;
    prec_bits += extra_guard;

    const unsigned int C_base = 640320u;
    mpz_class C3;
    mpz_ui_pow_ui(C3.get_mpz_t(), C_base, 3u);
    mpz_class C3_over_24 = C3 / 24;

#ifdef _OPENMP
    omp_set_num_threads(nthreads);
#endif

    std::vector<PQT> thread_results(nthreads);
    PQT result;

    // Manually split the work among threads using OpenMP 2.0
#pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        int num_threads = omp_get_num_threads();
        long long start = (thread_id * N) / num_threads;
        long long end = ((thread_id + 1) * N) / num_threads;

        if (start < end) {
            thread_results[thread_id] = binary_split(start, end, C3_over_24);
        }
    }

    // Combine the results from all threads sequentially in the main thread
    result = thread_results[0];
    for (int i = 1; i < nthreads; ++i) {
        if (thread_results[i].P != 0) { // Check if the thread had work to do
            combine(result, result, thread_results[i]);
        }
    }

    mpfr_set_default_prec((mpfr_prec_t)prec_bits);
    mpfr_t mp_Q, mp_T, mp_tmp, mp_num, mp_pi;
    mpfr_inits2((mpfr_prec_t)prec_bits, mp_Q, mp_T, mp_tmp, mp_num, mp_pi, (mpfr_ptr)0);

    mpfr_set_z(mp_Q, result.Q.get_mpz_t(), MPFR_RNDN);
    mpfr_set_z(mp_T, result.T.get_mpz_t(), MPFR_RNDN);

    mpfr_sqrt_ui(mp_tmp, 10005u, MPFR_RNDN);
    mpfr_mul_ui(mp_num, mp_Q, 426880u, MPFR_RNDN);
    mpfr_mul(mp_num, mp_num, mp_tmp, MPFR_RNDN);
    mpfr_div(mp_pi, mp_num, mp_T, MPFR_RNDN);

    char buf[64];
    mpfr_snprintf(buf, sizeof(buf), "%.15Rf", mp_pi);
    std::cout << "Pi (preview) = " << buf << std::endl;

    FILE* f = fopen(outfile, "w");
    if (!f) {
        std::cerr << "Cannot open output file: " << outfile << "\n";
        mpfr_clears(mp_Q, mp_T, mp_tmp, mp_num, mp_pi, (mpfr_ptr)0);
        return;
    }

    // Use mpfr_fprintf for robust, formatted output.
    // This prints pi with 'digits' number of digits after the decimal point.
    mpfr_fprintf(f, "%.*Rfd\n", (int)digits, mp_pi);

    fclose(f);

    mpfr_clears(mp_Q, mp_T, mp_tmp, mp_num, mp_pi, (mpfr_ptr)0);
}