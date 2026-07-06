#include "base_conv.hpp"
#include "bigint.hpp"
#include "ntt.hpp"
#include "timer.hpp"
#include "validator.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <gmp.h>
#include <iostream>
#include <omp.h>
#include <string>
#include <utility>
#include <vector>

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/resource.h>
#include <sys/time.h>
#endif

using namespace pi;

void log_event(const Timer &timer, const char *event,
               std::vector<std::pair<double, std::string>> &ev_list) {
  double t = timer.elapsed_seconds();
  printf("\n%.3f\t%s\n", t, event);
  fflush(stdout);
  ev_list.push_back({t, event});
}

int64_t parse_digits(std::string arg) {
  if (arg.empty())
    return 1000;
  char suffix = std::tolower(arg.back());
  int64_t multiplier = 1;
  bool has_suffix = true;
  if (suffix == 'k')
    multiplier = 1000;
  else if (suffix == 'm')
    multiplier = 1000000;
  else if (suffix == 'b')
    multiplier = 1000000000;
  else
    has_suffix = false;
  std::string num_part = has_suffix ? arg.substr(0, arg.size() - 1) : arg;
  try {
    return std::stoll(num_part) * multiplier;
  } catch (...) {
    return 1000;
  }
}

struct CPUMetrics {
  double user_time;
  double kernel_time;
};

CPUMetrics get_cpu_metrics() {
#ifdef _WIN32
  FILETIME creation_time, exit_time, kernel_time_ft, user_time_ft;
  if (GetProcessTimes(GetCurrentProcess(), &creation_time, &exit_time,
                      &kernel_time_ft, &user_time_ft)) {
    ULARGE_INTEGER kernel_ui, user_ui;
    kernel_ui.LowPart = kernel_time_ft.dwLowDateTime;
    kernel_ui.HighPart = kernel_time_ft.dwHighDateTime;
    user_ui.LowPart = user_time_ft.dwLowDateTime;
    user_ui.HighPart = user_time_ft.dwHighDateTime;
    return {(double)user_ui.QuadPart / 10000000.0,
            (double)kernel_ui.QuadPart / 10000000.0};
  }
#else
  struct rusage usage;
  if (getrusage(RUSAGE_SELF, &usage) == 0) {
    return {(double)usage.ru_utime.tv_sec +
                (double)usage.ru_utime.tv_usec / 1000000.0,
            (double)usage.ru_stime.tv_sec +
                (double)usage.ru_stime.tv_usec / 1000000.0};
  }
#endif
  return {0, 0};
}

int main(int argc, char *argv[]) {
  omp_set_max_active_levels(3);
  int64_t digits = 1000;
  if (argc > 1)
    digits = parse_digits(argv[1]);

  Timer total_timer;
  int64_t iterations = (double)digits / 14.1816459 + 10;
  std::vector<std::pair<double, std::string>> event_log;

  std::cout << "Program:               Pi-Calc v4"
            << std::endl;
  std::cout << "Algorithm:             Chudnovsky (1988)" << std::endl;
  std::cout << "Decimal Digits:        " << digits << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;

  std::cout << "Event Log:" << std::endl;
  log_event(total_timer, "Begin Computation", event_log);
  Timer comp_timer;

  BigInt P, Q, T;
  log_event(total_timer, "Step 1: Binary Splitting Start", event_log);
#pragma omp parallel
  {
#pragma omp single
    BigInt::binary_split(0, iterations, P, Q, T);
  }
  std::cout << std::endl;
  log_event(total_timer, "Step 1: Binary Splitting Finished", event_log);

  log_event(total_timer, "Step 2: Evaluation (Parallel)", event_log);
  int64_t guard = 64;
  mpz_t pi_z, num, sqrt_val, d10;
  mpz_init(pi_z);
  mpz_init(num);
  mpz_init(sqrt_val);
  mpz_init(d10);

  // High-performance evaluation
  BigInt::parallel_pow_ui(d10, 10, 2 * (digits + guard));
  mpz_mul_ui(d10, d10, 10005);
  mpz_sqrt(sqrt_val, d10);

  NTTMultiplier::multiply(num, Q.value, sqrt_val);
  mpz_mul_ui(num, num, 426880);
  mpz_tdiv_q(pi_z, num, T.value);

  mpz_t d10_guard;
  mpz_init(d10_guard);
  mpz_ui_pow_ui(d10_guard, 10, guard);
  mpz_tdiv_q(pi_z, pi_z, d10_guard);
  mpz_clear(d10_guard);
  log_event(total_timer, "Step 2: Evaluation Finished", event_log);

  double computation_time = comp_timer.elapsed_seconds();

  log_event(total_timer, "Step 3: Conversion & Writing Start", event_log);
  char *result_str = new char[digits + 5];
  BaseConverter::parallel_to_str(pi_z, digits + 1, result_str);

  FILE *f = fopen("pi.txt", "w");
  if (f) {
    fprintf(f, "%c.", result_str[0]);
    fwrite(result_str + 1, 1, digits, f);
    fclose(f);
  }
  log_event(total_timer, "Step 3: Conversion & Writing Finished", event_log);
  log_event(total_timer, "End Computation", event_log);

  double wall_time = total_timer.elapsed_seconds();
  CPUMetrics cpu = get_cpu_metrics();
  int threads = omp_get_max_threads();

  std::cout << "-----------------------------------------------" << std::endl;
  printf("Total Computation Time:    %.3f seconds\n", computation_time);
  printf("Start-to-End Wall Time:    %.3f seconds\n", wall_time);
  std::cout << std::endl;

  double user_util = (cpu.user_time / wall_time) * 100.0;
  double kernel_util = (cpu.kernel_time / wall_time) * 100.0;
  printf("CPU Utilization:           %.2f %%  +  %.2f %% kernel overhead\n",
         user_util, kernel_util);
  printf("Multi-core Efficiency:     %.2f %%  +  %.2f %% kernel overhead\n",
         user_util / threads, kernel_util / threads);
  std::cout << "-----------------------------------------------" << std::endl;

  // Run validation on generated decimal digits (after decimal point)
  ValidationResult val_res = PiValidator::validate(result_str + 1, digits);

  time_t now = time(nullptr);
  struct tm *tinfo = localtime(&now);
  char time_buf[64];
  strftime(time_buf, sizeof(time_buf), "%Y%m%d-%H%M%S", tinfo);

  std::string suffix_str;
  if (digits == 1000)
    suffix_str = "1K";
  else if (digits == 10000)
    suffix_str = "10K";
  else if (digits == 100000)
    suffix_str = "100K";
  else if (digits == 1000000)
    suffix_str = "1M";
  else if (digits == 10000000)
    suffix_str = "10M";
  else if (digits == 50000000)
    suffix_str = "50M";
  else if (digits == 100000000)
    suffix_str = "100M";
  else if (digits == 1000000000)
    suffix_str = "1B";
  else
    suffix_str = std::to_string(digits);

  std::string val_filename =
      "Pi - " + std::string(time_buf) + "(" + suffix_str + ").txt";

  PiValidator::write_validation_file(val_filename.c_str(), result_str + 1,
                                     digits, val_res, computation_time,
                                     wall_time, cpu.user_time, cpu.kernel_time,
                                     threads, event_log);

  delete[] result_str;
  mpz_clears(pi_z, num, sqrt_val, d10, NULL);
  return 0;
}
