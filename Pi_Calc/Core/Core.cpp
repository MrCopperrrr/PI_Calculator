#include "pch.h"
#include "Core.h"

// Temporarily disable warning C4146 just for the GMP header
#pragma warning(push)
#pragma warning(disable: 4146)
#include <gmp.h>
#pragma warning(pop)

#include <string>
#include <sstream>
#include <vector>
#include <thread>
#include <cmath>

// Constants for the Chudnovsky algorithm
const unsigned long A = 13591409;
const unsigned long B = 545140134;
const unsigned long C = 640320;
const unsigned long D = 12;

// Structure to hold the P, Q, T values for binary splitting
struct PQT {
    mpz_t p, q, t;
};

// This function computes a partial sum of the series for a given range [a, b)
// It is the core of the binary splitting algorithm.
void compute_bs(unsigned long a, unsigned long b, PQT& result) {
    mpz_init(result.p);
    mpz_init(result.q);
    mpz_init(result.t);

    if (b - a == 1) {
        // Base case of the recursion: compute for a single term 'a'
        if (a == 0) {
            mpz_set_ui(result.p, 1);
            mpz_set_ui(result.q, 1);
        }
        else {
            // P(a) = (6a-5)(2a-1)(6a-1)
            mpz_set_ui(result.p, 6 * a - 5);
            mpz_mul_ui(result.p, result.p, 2 * a - 1);
            mpz_mul_ui(result.p, result.p, 6 * a - 1);

            // Q(a) = C^3 * a^3 / 24
            mpz_t a_cubed;
            mpz_init(a_cubed);
            mpz_set_ui(a_cubed, a);
            mpz_pow_ui(a_cubed, a_cubed, 3);

            mpz_t C_cubed;
            mpz_init(C_cubed);
            mpz_set_ui(C_cubed, C);
            mpz_pow_ui(C_cubed, C_cubed, 3);

            mpz_mul(result.q, C_cubed, a_cubed);
            mpz_divexact_ui(result.q, result.q, 24);

            mpz_clear(a_cubed);
            mpz_clear(C_cubed);
        }

        // T(a) = (A + B*a) * P(a)
        mpz_set_ui(result.t, a);
        mpz_mul_ui(result.t, result.t, B);
        mpz_add_ui(result.t, result.t, A);
        mpz_mul(result.t, result.t, result.p);

        // T(a) *= (-1)^a
        if (a % 2 == 1) {
            mpz_neg(result.t, result.t);
        }
    }
    else {
        // Recursive step
        unsigned long mid = a + (b - a) / 2;
        PQT r1, r2;
        compute_bs(a, mid, r1);
        compute_bs(mid, b, r2);

        // Combine results from the two halves
        // P(a,b) = P(a,mid) * P(mid,b)
        // Q(a,b) = Q(a,mid) * Q(mid,b)
        // T(a,b) = T(a,mid)*Q(mid,b) + P(a,mid)*T(mid,b)
        mpz_mul(result.p, r1.p, r2.p);
        mpz_mul(result.q, r1.q, r2.q);

        mpz_t tmp1, tmp2;
        mpz_init(tmp1);
        mpz_init(tmp2);

        mpz_mul(tmp1, r1.t, r2.q);
        mpz_mul(tmp2, r1.p, r2.t);
        mpz_add(result.t, tmp1, tmp2);

        mpz_clear(r1.p); mpz_clear(r1.q); mpz_clear(r1.t);
        mpz_clear(r2.p); mpz_clear(r2.q); mpz_clear(r2.t);
        mpz_clear(tmp1); mpz_clear(tmp2);
    }
}


namespace PiCalculator
{
    std::string Chudnovsky::calculate(int digits, int num_threads)
    {
        // Set precision: digits * log2(10) + a small margin
        long precision_bits = (long)(digits * 3.3219280949) + 32;
        mpf_set_default_prec(precision_bits);

        // Determine the number of terms needed. Each term adds about 14.18 digits.
        unsigned long num_terms = (unsigned long)(digits / 14.18) + 2;

        std::vector<std::thread> threads;
        std::vector<PQT> results(num_threads);
        unsigned long terms_per_thread = num_terms / num_threads;

        for (int i = 0; i < num_threads; ++i) {
            unsigned long start = i * terms_per_thread;
            unsigned long end = (i == num_threads - 1) ? num_terms : start + terms_per_thread;
            if (start == end) continue; // Avoid launching threads that do no work
            threads.emplace_back(compute_bs, start, end, std::ref(results[i]));
        }

        for (auto& t : threads) {
            if (t.joinable()) t.join();
        }

        // Combine the results from all threads
        PQT final_res = results[0];
        for (size_t i = 1; i < threads.size(); ++i) {
            mpz_t p_temp, q_temp, t_temp;
            mpz_init_set(p_temp, final_res.p);
            mpz_init_set(q_temp, final_res.q);
            mpz_init_set(t_temp, final_res.t);

            mpz_mul(final_res.p, p_temp, results[i].p);
            mpz_mul(final_res.q, q_temp, results[i].q);

            mpz_t tmp1, tmp2;
            mpz_init(tmp1); mpz_init(tmp2);
            mpz_mul(tmp1, t_temp, results[i].q);
            mpz_mul(tmp2, p_temp, results[i].t);
            mpz_add(final_res.t, tmp1, tmp2);

            mpz_clear(p_temp); mpz_clear(q_temp); mpz_clear(t_temp);
            mpz_clear(tmp1); mpz_clear(tmp2);
            mpz_clear(results[i].p); mpz_clear(results[i].q); mpz_clear(results[i].t);
        }

        // Perform the final steps of the Chudnovsky formula
        // 1/pi = (12 * T) / (C^(3/2) * Q)
        mpf_t pi, C_sqrt, C_3_2, numerator, denominator;
        mpf_init(pi);
        mpf_init_set_ui(C_sqrt, C);
        mpf_init(C_3_2);
        mpf_init(numerator);
        mpf_init(denominator);

        mpf_sqrt(C_sqrt, C_sqrt);                      // C_sqrt = sqrt(C)
        mpf_mul_ui(C_3_2, C_sqrt, C);                  // C_3_2 = C * sqrt(C)

        mpf_set_z(numerator, final_res.t);
        mpf_mul_ui(numerator, numerator, D);           // numerator = 12 * T

        mpf_set_z(denominator, final_res.q);
        mpf_mul(denominator, denominator, C_3_2);      // denominator = Q * C^(3/2)

        mpf_div(pi, numerator, denominator);           // pi holds 1/pi for now
        mpf_ui_div(pi, 1, pi);                         // pi = 1 / (1/pi)

        // Format the output string
        mp_exp_t exp;
        char* str = mpf_get_str(NULL, &exp, 10, digits + 1, pi);
        std::string pi_str = "3." + std::string(str + 1);

        // Clean up
        free(str);
        mpz_clear(final_res.p);
        mpz_clear(final_res.q);
        mpz_clear(final_res.t);
        mpf_clear(pi);
        mpf_clear(C_sqrt);
        mpf_clear(C_3_2);
        mpf_clear(numerator);
        mpf_clear(denominator);

        return pi_str;
    }
}