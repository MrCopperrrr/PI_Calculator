#pragma once
#include "pqt.h"

// Ngưỡng để quyết định khi nào tạo task song song
static const long long TASK_THRESHOLD = 1024;

PQT binary_split_parallel(long long a, long long b, const mpz_class &C3_over_24) {
    // TRƯỜNG HỢP CƠ SỞ (LÁ CỦA CÂY ĐỆ QUY)
    if (b - a == 1) {
        PQT res;
        long long k = a;
        if (k == 0) {
            res.P = 1;
            res.Q = 1;
            res.T = 13591409;
            return res;
        }
        
        long long pk_ll = (6LL*k - 5) * (2LL*k - 1) * (6LL*k - 1);
        long long ak_ll = 13591409LL + 545140134LL * k;

        // SỬA LỖI 1: Tránh hàm tạo mơ hồ bằng cách dùng hàm C-style
        // mpz_set_si unambiguous, gán một số nguyên có dấu (signed integer).
        mpz_class k_mpz;
        mpz_set_si(k_mpz.get_mpz_t(), k);
        mpz_class qk = k_mpz * k_mpz * k_mpz * C3_over_24;

        __int128 tk_128 = (__int128)pk_ll * ak_ll;
        if (k & 1) {
            tk_128 = -tk_128;
        }

        // SỬA LỖI 2: Tránh hàm tạo mơ hồ tương tự như trên
        mpz_set_si(res.P.get_mpz_t(), pk_ll);
        res.Q = qk;

        bool is_negative = tk_128 < 0;
        unsigned __int128 val_abs = is_negative ? -tk_128 : tk_128;

        unsigned long long limbs[2];
        limbs[0] = (unsigned long long)(val_abs >> 64);
        limbs[1] = (unsigned long long)val_abs;

        // SỬA LỖI 3: Dùng kiểu con trỏ mpz_ptr thay vì kiểu mảng mpz_t
        // get_mpz_t() trả về một con trỏ, không phải một mảng.
        mpz_ptr t_mpz = res.T.get_mpz_t();

        mpz_import(t_mpz, 2, 1, sizeof(unsigned long long), 0, 0, limbs);

        if (is_negative) {
            mpz_neg(t_mpz, t_mpz);
        }
        
        return res;
    }

    // PHẦN ĐỆ QUY SONG SONG (giữ nguyên, đã đúng)
    long long m = (a + b) >> 1;
    PQT res, L, R;

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