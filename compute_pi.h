#include "binary_split.h"
#include <mpfr.h>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>

using namespace std;

void compute_pi_bs_parallel(long long digits, int nthreads, const string &outfile) {
    if (digits < 1) digits = 1;
    if (nthreads < 1) nthreads = 1;

    long long N = (digits + 13) / 14;

    long long prec_bits = (long long)ceil(digits * 3.3219280948873626);
    long long extra_guard = 512 + (long long)std::min( (long long)4096, (long long)(digits / 1000000) );
    prec_bits += extra_guard;

    // C^3 và C^3/24
    const unsigned int C_base = 640320u;
    mpz_class C3;
    mpz_ui_pow_ui(C3.get_mpz_t(), C_base, 3u); // C^3
    mpz_class C3_over_24 = C3 / mpz_class(24);

    #ifdef _OPENMP
    omp_set_num_threads(nthreads);
    #endif

    PQT result;
    #pragma omp parallel
    {
        #pragma omp single
        {
            result = binary_split_parallel(0, N, C3_over_24);
        }
    }

    // MPFR: pi = (Q * 426880 * sqrt(10005)) / T
    mpfr_set_default_prec((mpfr_prec_t)prec_bits);
    mpfr_t mp_Q, mp_T, mp_tmp, mp_num, mp_pi;
    mpfr_init2(mp_Q,   prec_bits);
    mpfr_init2(mp_T,   prec_bits);
    mpfr_init2(mp_tmp, prec_bits);
    mpfr_init2(mp_num, prec_bits);
    mpfr_init2(mp_pi,  prec_bits);

    mpfr_set_z(mp_Q, result.Q.get_mpz_t(), MPFR_RNDN);
    mpfr_set_z(mp_T, result.T.get_mpz_t(), MPFR_RNDN);

    mpfr_set_ui(mp_tmp, 10005u, MPFR_RNDN);
    mpfr_sqrt(mp_tmp, mp_tmp, MPFR_RNDN);             // sqrt(10005)

    mpfr_set_z(mp_num, result.Q.get_mpz_t(), MPFR_RNDN);
    mpfr_mul_ui(mp_num, mp_num, 426880u, MPFR_RNDN);  // Q * 426880
    mpfr_mul(mp_num, mp_num, mp_tmp, MPFR_RNDN);      // * sqrt(10005)

    mpfr_div(mp_pi, mp_num, mp_T, MPFR_RNDN);         // pi

    {
        mpfr_t pi_debug;
        mpfr_init2(pi_debug, prec_bits);
        mpfr_set(pi_debug, mp_pi, MPFR_RNDN);

        char buf[64];
        // In ra 15 chữ số để chắc chắn
        mpfr_snprintf(buf, sizeof(buf), "%.15Rf", pi_debug);
        cout << "Pi (preview) = " << buf << endl;

        mpfr_clear(pi_debug);
    }

    // ===== Streaming output =====
    FILE *f = fopen(outfile.c_str(), "w");
    if (!f) {
        cerr << "Cannot open output file: " << outfile << "\n";
        mpfr_clears(mp_Q, mp_T, mp_tmp, mp_num, mp_pi, (mpfr_ptr)0);
        return;
    }

    mpz_class zint;
    mpfr_get_z(zint.get_mpz_t(), mp_pi, MPFR_RNDZ);   // truncate toward 0
    std::string int_str = zint.get_str(10);

    fputs(int_str.c_str(), f);
    fputc('.', f);

    // frac = pi - floor(pi)
    mpfr_t frac, tmp;
    mpfr_init2(frac, prec_bits);
    mpfr_init2(tmp,  prec_bits);

    mpfr_set_z(tmp, zint.get_mpz_t(), MPFR_RNDN);
    mpfr_sub(frac, mp_pi, tmp, MPFR_RNDN);            // 0 <= frac < 1

    const int BLOCK = 5000000; // 500k digits per block
    mpfr_t pow10_10k, pow10_var;
    mpfr_init2(pow10_10k, prec_bits);
    mpfr_init2(pow10_var,  prec_bits);

    mpz_t pow10_block_z;
    mpz_init(pow10_block_z);
    if (BLOCK > 0) mpz_ui_pow_ui(pow10_block_z, 10u, (unsigned long)BLOCK);
    mpfr_set_z(pow10_10k, pow10_block_z, MPFR_RNDN);

    mpz_t zblock;
    mpz_init(zblock);

    long long remaining = digits; 
    while (remaining > 0) {
        int take = (remaining >= BLOCK) ? BLOCK : (int)remaining;

        if (take == BLOCK) {
            mpfr_mul(tmp, frac, pow10_10k, MPFR_RNDN);
        } else {
            mpfr_set_ui(pow10_var, 10u, MPFR_RNDN);
            mpfr_pow_ui(pow10_var, pow10_var, (unsigned long)take, MPFR_RNDN);
            mpfr_mul(tmp, frac, pow10_var, MPFR_RNDN);
        }

        mpfr_get_z(zblock, tmp, MPFR_RNDZ); // truncate

        std::string blk = mpz_class(zblock).get_str(10);
        if ((int)blk.size() < take) {
            std::string pad(take - (int)blk.size(), '0');
            fwrite(pad.data(), 1, pad.size(), f);
        }
        fwrite(blk.data(), 1, blk.size(), f);

        mpfr_set_z(mp_tmp, zblock, MPFR_RNDN);
        mpfr_sub(frac, tmp, mp_tmp, MPFR_RNDN);

        remaining -= take;

        fflush(f);
    }

    fputc('\n', f);
    fclose(f);

    mpz_clear(zblock);
    mpfr_clears(mp_Q, mp_T, mp_tmp, mp_num, mp_pi, frac, tmp, pow10_10k, pow10_var, (mpfr_ptr)0);
}
