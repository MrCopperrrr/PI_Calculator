#include "pch.h"
#include "Core.h"
#include <gmp.h>
#include <string>
#include <sstream>
#include <vector>
#include <thread>
#include <cmath>
#include <mutex>

// Constants for the Chudnovsky algorithm
#define A 13591409
#define B 545140134
#define C 640320
#define D 12

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
        // Base case of the recursion
        mpz_set_ui(result.p, (2 * b - 1));
        mpz_mul_ui(result.p, result.p, (6 * b - 1));
        mpz_mul_ui(result.p, result.p, (6 * b - 5));

        mpz_set_ui(result.q, C);
        mpz_mul_ui(result.q, result.q, C);
        mpz_mul_ui(result.q, result.q, C);
        mpz_divexact_ui(result.q, result.q, 24);

        mpz_set_si(result.t, (A + (long long)B * b));
        if (b % 2 == 1) {
            mpz_neg(result.t, result.t);
        }
        mpz_mul(result.t, result.t, result.p);
    }
    else {
        // Recursive step
        unsigned long mid = a + (b - a) / 2;
        PQT r1, r2;
        compute_bs(a, mid, r1);
        compute_bs(mid, b, r2);

        // Combine results from the two halves
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
        long precision_bits = (long)(digits * 3.3219280949) + 16;
        mpf_set_default_prec(precision_bits);

        // Determine the number of terms needed for the series.
        // Each term adds about 14.18 digits of precision.
        unsigned long num_terms = (unsigned long)(digits / 14.18) + 2;

        std::vector<std::thread> threads;
        std::vector<PQT> results(num_threads);
        unsigned long terms_per_thread = num_terms / num_threads;

        // Launch threads to compute parts of the series in parallel
        for (int i = 0; i < num_threads; ++i) {
            unsigned long start = i * terms_per_thread;
            unsigned long end = (i == num_threads - 1) ? num_terms : start + terms_per_thread;
            threads.emplace_back(compute_bs, start, end, std::ref(results[i]));
        }

        // Wait for all threads to finish
        for (auto& t : threads) {
            t.join();
        }

        // --- Combine the results from all threads ---
        PQT final_res = results[0];
        for (int i = 1; i < num_threads; ++i) {
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

        // --- Perform the final steps of the Chudnovsky formula ---
        mpf_t pi, temp_q, temp_t, C_sqrt;
        mpf_init(pi);
        mpf_init(temp_q);
        mpf_init(temp_t);
        mpf_init_set_ui(C_sqrt, C);

        mpf_sqrt(C_sqrt, C_sqrt); // sqrt(C)

        mpf_set_z(temp_q, final_res.q); // Q value
        mpf_set_z(temp_t, final_res.t); // T value

        // pi = (Q * sqrt(C)) / (A*Q + T)
        mpf_mul_ui(pi, temp_q, A);
        mpf_add(pi, pi, temp_t);
        mpf_ui_div(pi, D, pi); // This is part of the formula, D*...

        mpf_mul(pi, pi, temp_q);
        mpf_mul(pi, pi, C_sqrt);

        mpf_t inv_pi;
        mpf_init(inv_pi);
        mpf_ui_div(inv_pi, 1, pi);

        // --- Format the output string ---
        // Request digits + 2 characters ('3' and the null terminator)
        // mpf_get_str allocates memory, which must be freed.
        mp_exp_t exp;
        char* str = mpf_get_str(NULL, &exp, 10, digits + 1, inv_pi);
        std::string pi_str = "3." + std::string(str + 1);

        // Clean up all GMP variables
        free(str);
        mpz_clear(final_res.p);
        mpz_clear(final_res.q);
        mpz_clear(final_res.t);
        mpf_clear(pi);
        mpf_clear(temp_q);
        mpf_clear(temp_t);
        mpf_clear(C_sqrt);
        mpf_clear(inv_pi);

        return pi_str;
    }
}