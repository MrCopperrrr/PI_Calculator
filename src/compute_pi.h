#pragma once
#include "binary_split.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>
#include <thread>      // Cần cho std::thread
#include <queue>       // Cần cho std::queue
#include <mutex>       // Cần cho std::mutex
#include <condition_variable> // Cần cho std::condition_variable

#include <gmpxx.h>

// =========================================================================
// LỚP HÀNG ĐỢI AN TOÀN LUỒNG (THREAD-SAFE QUEUE)
// =========================================================================
template<typename T>
class ThreadSafeQueue {
public:
    void push(T value) {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_queue.push(std::move(value));
        m_cond.notify_one();
    }

    bool pop(T& value) {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_cond.wait(lock, [this]{ return !m_queue.empty() || m_finished; });
        if (m_queue.empty() && m_finished) {
            return false; // Hàng đợi trống và đã kết thúc
        }
        value = std::move(m_queue.front());
        m_queue.pop();
        return true;
    }

    void finish() {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_finished = true;
        m_cond.notify_all();
    }

private:
    std::queue<T> m_queue;
    std::mutex m_mutex;
    std::condition_variable m_cond;
    bool m_finished = false;
};

// Kiểu dữ liệu để truyền giữa các luồng
struct PiBlock {
    std::string data;
    int expected_length;
};

// =========================================================================
// HÀM CONSUMER (GHI FILE)
// =========================================================================
void file_writer_task(ThreadSafeQueue<PiBlock>& queue, const std::string& outfile) {
    FILE *f = fopen(outfile.c_str(), "a"); // Mở file ở chế độ "append"
    if (!f) {
        std::cerr << "Writer thread: Cannot open output file.\n";
        return;
    }

    PiBlock block;
    while (queue.pop(block)) {
        if ((int)block.data.length() < block.expected_length) {
            std::string padding(block.expected_length - block.data.length(), '0');
            fputs(padding.c_str(), f);
        }
        fputs(block.data.c_str(), f);
        fflush(f);
    }
    
    fputc('\n', f);
    fclose(f);
}

// =========================================================================
// HÀM TÍNH PI CHÍNH (ĐÃ SỬA ĐỔI)
// =========================================================================
void compute_pi_bs_parallel(long long digits, int nthreads, const std::string &outfile) {
    if (digits < 1) digits = 1;
    if (nthreads < 1) nthreads = 1;

    // ... (Phần tính toán P, Q, T không đổi) ...
    long long N = (digits / 14.18164746) + 1;
    std::cout << "Calculating " << N << " terms for " << digits << " digits.\n";
    const unsigned int C_base = 640320u;
    mpz_class C3;
    mpz_ui_pow_ui(C3.get_mpz_t(), C_base, 3u);
    mpz_class C3_over_24 = C3 / 24;
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

    // ... (Phần tính Pi bằng mpf_class không đổi) ...
    long long final_prec_bits = (long long)ceil(digits * 3.3219280948873626) + 64;
    std::cout << "Required precision: " << final_prec_bits << " bits.\n";
    mpf_set_default_prec(final_prec_bits);
    mpf_class num;
    mpf_class const_sqrt;
    mpf_sqrt_ui(const_sqrt.get_mpf_t(), 10005);
    num = mpf_class(result.Q);
    num *= 426880;
    num *= const_sqrt;
    mpf_class den(result.T);
    mpf_class pi = num / den;

    // === BẮT ĐẦU PHẦN PIPELINING ===

    // 1. Chuẩn bị file output (ghi phần nguyên trước)
    {
        FILE *f = fopen(outfile.c_str(), "w");
        if (!f) {
            std::cerr << "Cannot open output file: " << outfile << "\n";
            return;
        }
        mpz_class zint(pi);
        std::string int_str = zint.get_str();
        fputs(int_str.c_str(), f);
        fputc('.', f);
        fclose(f);
    }

    // 2. Khởi tạo hàng đợi và luồng ghi file (Consumer)
    ThreadSafeQueue<PiBlock> writer_queue;
    std::thread writer_thread(file_writer_task, std::ref(writer_queue), std::ref(outfile));

    // 3. Luồng chính trở thành Producer
    mpf_class frac = pi - mpz_class(pi);
    // ... (Logic chọn block_size linh hoạt không đổi) ...
    // 1. Đặt ra giới hạn hợp lý cho block size (giữ nguyên)
    const int MIN_BLOCK_SIZE = 1000000;      // Một triệu
    const int MAX_BLOCK_SIZE = 25000000;     // Hai mươi lăm triệu

    // 2. Xây dựng heuristic thông minh hơn dựa trên cả `digits` và `nthreads`
    const double BASE_THREADS = 16.0; // Chọn 16 luồng làm cơ sở
    // Tính hệ số điều chỉnh, giới hạn nó trong một khoảng hợp lý (ví dụ: 0.5 đến 2.0)
    // để tránh các giá trị quá cực đoan khi nthreads rất thấp hoặc rất cao.
    double adjustment_factor = std::sqrt(BASE_THREADS / nthreads);
    adjustment_factor = std::max(0.5, std::min(adjustment_factor, 2.0));

    // Số lần lặp mục tiêu, được điều chỉnh bởi hệ số
    double target_iterations = 100.0 * adjustment_factor;

    // Tính block_size lý tưởng
    int block_size = digits / target_iterations;

    // 3. Giới hạn block_size trong khoảng [MIN, MAX] (giữ nguyên)
    block_size = std::max(MIN_BLOCK_SIZE, std::min(block_size, MAX_BLOCK_SIZE));
    
    // Nếu tổng số chữ số nhỏ hơn cả MIN_BLOCK_SIZE, thì chỉ dùng block_size = digits
    if (digits < MIN_BLOCK_SIZE) {
        block_size = digits;
    }

    std::cout << "Using dynamic block size (adjusted for " << nthreads << " threads): " 
              << block_size << " digits per block.\n";

    // ... (Logic tính pow10_block không đổi) ...
    std::cout << "Pre-calculating 10^" << block_size << " for block processing...\n";
    mpz_class pow10_int;
    mpz_ui_pow_ui(pow10_int.get_mpz_t(), 10, block_size);
    mpf_class pow10_block(pow10_int);
    std::cout << "Pre-calculation complete.\n";

    long long remaining = digits;
    while (remaining > 0) {
        int take = (remaining >= block_size) ? block_size : (int)remaining;
        
        // --- PRODUCER: TÍNH TOÁN ---
        if (take < block_size) {
            mpz_class pow10_temp_int;
            mpz_ui_pow_ui(pow10_temp_int.get_mpz_t(), 10, take);
            frac *= mpf_class(pow10_temp_int);
        } else {
            frac *= pow10_block;
        }
        mpz_class zblock(frac);
        
        // --- PRODUCER: ĐẨY VÀO HÀNG ĐỢI ---
        PiBlock block_to_write;
        block_to_write.data = zblock.get_str();
        block_to_write.expected_length = take;
        writer_queue.push(block_to_write);

        frac -= zblock;
        remaining -= take;
    }

    // 4. Báo cho luồng ghi file biết là đã hết việc
    writer_queue.finish();

    // 5. Chờ luồng ghi file hoàn thành công việc của nó
    writer_thread.join();
}