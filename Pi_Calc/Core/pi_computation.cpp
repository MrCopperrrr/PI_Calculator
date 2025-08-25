#include "core_api.h"
#include "binary_split.h"
#include <mpfr.h>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

// This is the SERIAL implementation of the binary split algorithm.
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

        // ====================================================================
        // FINAL, FINAL ROBUST IMPLEMENTATION
        // The fix is here: explicitly cast the long long 'k' to a 'long'
        // to resolve the constructor ambiguity.
        // ====================================================================
        mpz_class K(static_cast<long>(k));

        // p(k) = (6k-5)(2k-1)(6k-1)
        mpz_class pk = (mpz_class(6) * K - 5);
        pk *= (mpz_class(2) * K - 1);
        pk *= (mpz_class(6) * K - 1);

        // q(k) = k^3 * (C^3 / 24)
        mpz_class qk = K * K * K;
        qk *= C3_over_24;

        // a(k) = 13591409 + 545140134 * k
        mpz_class ak = mpz_class(13591409) + mpz_class(545140134) * K;

        mpz_class tk = pk * ak;
        if (k & 1) {
            tk = -tk;
        }

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
    if (N == 0) N = 1;

    const unsigned int C_base = 640320u;
    mpz_class C3;
    mpz_ui_pow_ui(C3.get_mpz_t(), C_base, 3u);
    mpz_class C3_over_24 = C3 / 24;

#ifdef _OPENMP
    omp_set_num_threads(nthreads);
#endif

    std::vector<PQT> thread_results(nthreads);
    PQT result;

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

    // Robustly combine thread results
    int first_valid_result_idx = -1;
    for (int i = 0; i < nthreads; ++i) {
        if (thread_results[i].P != 0) {
            first_valid_result_idx = i;
            break;
        }
    }
    if (first_valid_result_idx != -1) {
        result = thread_results[first_valid_result_idx];
        for (int i = first_valid_result_idx + 1; i < nthreads; ++i) {
            if (thread_results[i].P != 0) {
                combine(result, result, thread_results[i]);
            }
        }
    }
    else {
        result = binary_split(0, 1, C3_over_24);
    }

    // ====================================================================
    // PRECISION FIX: The precision must be large enough for the
    // intermediate integers AND the final result.
    // ====================================================================
    long long final_prec_bits = (long long)ceil(digits * 3.3219280948873626);

    // Calculate bits needed for the huge integers Q and T
    long long integer_prec_bits = std::max(
        mpz_sizeinbase(result.Q.get_mpz_t(), 2),
        mpz_sizeinbase(result.T.get_mpz_t(), 2)
    );

    // Use the larger of the two precisions, plus a safety margin
    long long prec_bits = std::max(final_prec_bits, integer_prec_bits) + 128;
    // ====================================================================
    // END OF FIX
    // ====================================================================

    mpfr_set_default_prec((mpfr_prec_t)prec_bits);
    mpfr_t mp_Q, mp_T, mp_tmp, mp_num, mp_pi;
    mpfr_inits2((mpfr_prec_t)prec_bits, mp_Q, mp_T, mp_tmp, mp_num, mp_pi, (mpfr_ptr)0);

    // These conversions are now safe because prec_bits is large enough
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

    mpfr_fprintf(f, "%.*Rfd\n", (int)digits, mp_pi);

    fclose(f);

    mpfr_clears(mp_Q, mp_T, mp_tmp, mp_num, mp_pi, (mpfr_ptr)0);
}