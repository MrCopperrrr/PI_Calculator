#include "pch.h"
#include "Core.h"
#include "..\GPU\gpu_multiplier.h" // The header for your GPU functions

#include <gmp.h>
#include <string>
#include <sstream>
#include <vector>
#include <thread>
#include <cmath>
#include <algorithm> // For std::max
#include <iostream>  // <--- FIX #1: Added this header for std::cout

// --- HYBRID CPU+GPU CONFIGURATION ---
constexpr size_t GPU_MULTIPLY_THRESHOLD_BITS = 65536;

// --- GMP <-> RAW ARRAY CONVERSION HELPERS ---
void mpz_to_uint_vector(const mpz_t gmp_num, std::vector<unsigned int>& vec) {
    if (mpz_sgn(gmp_num) == 0) { vec.clear(); return; }
    size_t num_bytes = (mpz_sizeinbase(gmp_num, 2) + 7) / 8;
    size_t vec_size = (num_bytes + sizeof(unsigned int) - 1) / sizeof(unsigned int);
    vec.resize(vec_size);
    size_t count;
    mpz_export(vec.data(), &count, -1, sizeof(unsigned int), 0, 0, gmp_num);
    vec.resize(count);
}

void uint_vector_to_mpz(const std::vector<unsigned int>& vec, mpz_t gmp_num) {
    if (vec.empty()) { mpz_set_ui(gmp_num, 0); return; }
    mpz_import(gmp_num, vec.size(), -1, sizeof(unsigned int), 0, 0, vec.data());
}

// --- HYBRID MULTIPLICATION FUNCTION ---
void hybrid_mul(mpz_t result, const mpz_t op1, const mpz_t op2) {
    size_t op1_bits = mpz_sizeinbase(op1, 2);
    size_t op2_bits = mpz_sizeinbase(op2, 2);
    if (std::max(op1_bits, op2_bits) > GPU_MULTIPLY_THRESHOLD_BITS) {
        std::vector<unsigned int> op1_vec, op2_vec;
        mpz_to_uint_vector(op1, op1_vec);
        mpz_to_uint_vector(op2, op2_vec);
        int result_size = static_cast<int>(op1_vec.size() + op2_vec.size() + 1);
        std::vector<unsigned int> result_vec(result_size);
        multiply_on_gpu(op1_vec.data(), static_cast<int>(op1_vec.size()),
            op2_vec.data(), static_cast<int>(op2_vec.size()),
            result_vec.data(), result_size);
        result_vec.resize(result_size);
        uint_vector_to_mpz(result_vec, result);
    }
    else {
        mpz_mul(result, op1, op2);
    }
}

// --- CHUDNOVSKY ALGORITHM IMPLEMENTATION ---
#define A 13591409
#define B 545140134
#define C 640320
#define D 12
struct PQT { mpz_t p, q, t; };

void compute_bs(unsigned long a, unsigned long b, PQT& result) {
    mpz_init(result.p); mpz_init(result.q); mpz_init(result.t);
    if (b - a == 1) {
        mpz_set_ui(result.p, (2 * b - 1));
        mpz_mul_ui(result.p, result.p, (6 * b - 1));
        mpz_mul_ui(result.p, result.p, (6 * b - 5));
        mpz_set_ui(result.q, C);
        mpz_pow_ui(result.q, result.q, 3);
        mpz_divexact_ui(result.q, result.q, 24);
        mpz_set_si(result.t, (A + (long long)B * b));
        if (b % 2 == 1) { mpz_neg(result.t, result.t); }
        mpz_mul(result.t, result.t, result.p);
    }
    else {
        unsigned long mid = a + (b - a) / 2;
        PQT r1, r2;
        compute_bs(a, mid, r1);
        compute_bs(mid, b, r2);
        hybrid_mul(result.p, r1.p, r2.p);
        hybrid_mul(result.q, r1.q, r2.q);
        mpz_t tmp1, tmp2;
        mpz_init(tmp1); mpz_init(tmp2);
        hybrid_mul(tmp1, r1.t, r2.q);
        hybrid_mul(tmp2, r1.p, r2.t);
        mpz_add(result.t, tmp1, tmp2);
        mpz_clear(r1.p); mpz_clear(r1.q); mpz_clear(r1.t);
        mpz_clear(r2.p); mpz_clear(r2.q); mpz_clear(r2.t);
        mpz_clear(tmp1); mpz_clear(tmp2);
    }
}

namespace PiCalculator {
    std::string Chudnovsky::calculate(int digits, int num_threads) {
        long precision_bits = (long)(digits * 3.3219280949) + 16;
        mpf_set_default_prec(precision_bits);
        unsigned long num_terms = (unsigned long)(digits / 14.18) + 2;
        int effective_num_threads = std::max(1, std::min((int)num_terms, num_threads));
        if (num_threads != effective_num_threads) {
            std::cout << "[INFO] Using " << effective_num_threads << " threads (capped from " << num_threads << ") for efficiency." << std::endl;
        }
        std::vector<std::thread> threads;
        std::vector<PQT> results(effective_num_threads);
        unsigned long terms_per_thread = num_terms / effective_num_threads;
        for (int i = 0; i < effective_num_threads; ++i) {
            unsigned long start = i * terms_per_thread;
            unsigned long end = (i == effective_num_threads - 1) ? num_terms : start + terms_per_thread;
            if (start >= end) continue;
            threads.emplace_back(compute_bs, start, end, std::ref(results[i]));
        }
        for (auto& t : threads) { t.join(); }
        if (results.empty() || threads.empty()) {
            return "Error: No work was performed. Check input digits.";
        }
        PQT final_res = results[0];
        for (size_t i = 1; i < results.size(); ++i) {
            mpz_t p_temp, q_temp, t_temp;
            mpz_init_set(p_temp, final_res.p);
            mpz_init_set(q_temp, final_res.q);
            mpz_init_set(t_temp, final_res.t);
            hybrid_mul(final_res.p, p_temp, results[i].p);
            hybrid_mul(final_res.q, q_temp, results[i].q);
            mpz_t tmp1, tmp2;
            mpz_init(tmp1); mpz_init(tmp2);
            hybrid_mul(tmp1, t_temp, results[i].q);
            hybrid_mul(tmp2, p_temp, results[i].t);
            mpz_add(final_res.t, tmp1, tmp2);
            mpz_clear(p_temp); mpz_clear(q_temp); mpz_clear(t_temp);
            mpz_clear(tmp1); mpz_clear(tmp2);
            mpz_clear(results[i].p); mpz_clear(results[i].q); mpz_clear(results[i].t);
        }
        mpf_t pi, C_sqrt;
        mpf_init(pi); mpf_init_set_ui(C_sqrt, C);
        mpf_sqrt(C_sqrt, C_sqrt);
        mpf_t final_q, final_t;
        mpf_init(final_q); mpf_init(final_t);

        // --- FIX #2: Changed mpz_get_f to the correct mpf_set_z ---
        mpf_set_z(final_q, final_res.q);
        mpf_set_z(final_t, final_res.t);

        mpf_mul(pi, final_q, C_sqrt);
        mpf_t temp_sum;
        mpf_init(temp_sum);
        mpf_mul_ui(temp_sum, final_q, A);
        mpf_add(temp_sum, temp_sum, final_t);
        mpf_div(pi, pi, temp_sum);
        mpf_mul_ui(pi, pi, D);
        mp_exp_t exp;
        char* str = mpf_get_str(NULL, &exp, 10, digits + 1, pi);
        std::string pi_str = "3." + std::string(str + 1);
        free(str);
        mpz_clear(final_res.p); mpz_clear(final_res.q); mpz_clear(final_res.t);
        mpf_clear(pi); mpf_clear(C_sqrt); mpf_clear(final_q); mpf_clear(final_t); mpf_clear(temp_sum);
        return pi_str;
    }
}