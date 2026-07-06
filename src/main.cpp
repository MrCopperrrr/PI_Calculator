#include "base_conv.hpp"
#include "bigint.hpp"
#include "ntt.hpp"
#include "timer.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <gmp.h>
#include <iostream>
#include <omp.h>
#include <string>
#include <vector>


#ifdef _WIN32
#include <windows.h>
#else
#include <sys/resource.h>
#include <sys/time.h>
#endif

using namespace pi;

#include <ctime>
#include <iomanip>
#include <fstream>
#include <sstream>

using namespace pi;

void log_event(const Timer &timer, const char *event) {
  printf("\n%.3f\t%s\n", timer.elapsed_seconds(), event);
  fflush(stdout);
}

// Helper to get formatted current time
std::string get_timestamp() {
  std::time_t now = std::time(nullptr);
  char buf[20];
  std::strftime(buf, sizeof(buf), "%Y%m%d-%H%M%S", std::localtime(&now));
  return std::string(buf);
}

void write_validation_report(int64_t digits, double comp_time, double wall_time,
                             double user_util, double kernel_util,
                             const std::vector<std::pair<std::string, double>>& events,
                             const char* result_str, size_t actual_len) {
  time_t now = time(0);
  tm* ltm = localtime(&now);
  char timestamp[20];
  strftime(timestamp, 20, "%Y%m%d-%H%M%S", ltm);
  std::string filename = "Validation - Pi - " + std::string(timestamp) + ".txt";

  std::ofstream f(filename);
  if (!f.is_open()) return;

  f << "Pi-Calc Validation File - High Precision Computation Result\n";
  f << "-----------------------------------------------------------\n\n";
  f << "Program:               Pi-Calc v2.5 (Extreme Edition)\n";
  f << "Algorithm:             Chudnovsky (1988)\n";
  f << "Decimal Digits:        " << digits << "\n";
  f << "Threading Mode:        OpenMP Tasking -> " << omp_get_max_threads() << " threads\n\n";
  f << "Computation Metrics:\n";
  f << "Total Computation Time:    " << std::fixed << std::setprecision(3) << comp_time << " seconds\n";
  f << "Start-to-End Wall Time:    " << wall_time << " seconds\n";
  f << "CPU Utilization:           " << user_util << " %  +  " << kernel_util << " % kernel overhead\n\n";
  f << "Final String Length:       " << actual_len << " (Expected " << digits + 1 << ")\n\n";

  f << "Event Log:\n";
  for (const auto& ev : events) {
    f << std::setw(8) << std::left << ev.second << ev.first << "\n";
  }
  f << "\n";

  f << "Decimal Digits Samples:\n";
  // Sample at point s: show 50 digits ending at position s
  std::vector<int64_t> samples = {50, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000};
  for (int64_t s : samples) {
    if (s <= digits && (size_t)s < actual_len) {
        f << std::setw(12) << std::left << (s == 50 ? "3. 51" : std::to_string(s)) << ":  ";
        for (int i = 0; i < 50; ++i) {
            if (i > 0 && i % 10 == 0) f << " ";
            int64_t idx = s - 50 + 1 + i; 
            if (idx < (int64_t)actual_len) f << result_str[idx];
            else f << "?";
        }
        f << "\n";
    }
  }
  
  if (actual_len > 50) {
      f << std::setw(12) << std::left << "End" << ":  ";
      for (int i = 0; i < 50; ++i) {
          if (i > 0 && i % 10 == 0) f << " ";
          int64_t idx = actual_len - 51 + i;
          if (idx >= 0) f << result_str[idx];
      }
      f << "\n";
  }

  f << "\n-----------------------------------------------------------\n";
  f << "Validation Status: Spot Check Recommended with y-cruncher logs.\n";
  f.close();
  
  std::cout << "Validation report saved to: " << filename << std::endl;
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
  std::vector<std::pair<std::string, double>> event_history;
  auto record_event = [&](const char* name) {
      double t = total_timer.elapsed_seconds();
#pragma omp critical(event_log)
      {
          event_history.push_back({name, t});
          log_event(total_timer, name);
      }
  };

  int64_t iterations = (double)digits / 14.181647 + 200;
  int64_t guard = 256;

  std::cout << "Program:               Pi-Calc (Version 3.0)"
            << std::endl;
  std::cout << "Algorithm:             Chudnovsky (1988)" << std::endl;
  std::cout << "Decimal Digits:        " << digits << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;

  std::cout << "Event Log:" << std::endl;
  record_event("Begin Computation");
  Timer comp_timer;

  BigInt P, Q, T;
  record_event("Step 1: Binary Splitting Start");
#pragma omp parallel
  {
#pragma omp single
    BigInt::binary_split(0, iterations, P, Q, T);
  }
  std::cout << std::endl;
  record_event("Step 1: Binary Splitting Finished");
  
  NTTMultiplier::use_hybrid = true; // Switch to multi-core strategy for Step 2
  record_event("Step 2: Evaluation (Parallel)");
  
  mpz_t pi_z, num, sqrt_val, d10;
  mpz_init(pi_z);
  mpz_init(num);
  mpz_init(sqrt_val);
  mpz_init(d10);

  record_event("Step 2.1: Power of 10 Start");
  mpz_ui_pow_ui(d10, 10, 2 * (digits + guard));
  mpz_mul_ui(d10, d10, 10005);
  record_event("Step 2.1: Power of 10 Finished");

  record_event("Step 2.2: Square Root Start");
  BigInt::parallel_sqrt(sqrt_val, d10);
  record_event("Step 2.2: Square Root Finished");

  record_event("Step 2.3: Multiplier Start");
  NTTMultiplier::multiply(num, Q.value, sqrt_val);
  mpz_mul_ui(num, num, 426880);
  record_event("Step 2.3: Multiplier Finished");

  record_event("Step 2.4: Final Division Start");
  BigInt::parallel_div(pi_z, num, T.value);
  
  mpz_ui_pow_ui(d10, 10, guard);
  mpz_tdiv_q(pi_z, pi_z, d10);
  record_event("Step 2.4: Final Division Finished");

  record_event("Step 2: Evaluation Finished");

  double computation_time = comp_timer.elapsed_seconds();

  record_event("Step 3: Conversion & Writing Start");
  char *result_str = new char[digits + 5];
  BaseConverter::parallel_to_str(pi_z, digits + 1, result_str);
  size_t actual_len = digits + 1;

  FILE *f = fopen("pi.txt", "w");
  if (f) {
    fprintf(f, "%c.", result_str[0]);
    fwrite(result_str + 1, 1, digits, f);
    fclose(f);
  }
  record_event("Step 3: Conversion & Writing Finished");
  record_event("End Computation");

  double wall_time = total_timer.elapsed_seconds();
  CPUMetrics cpu = get_cpu_metrics();
  int threads = omp_get_max_threads();

  double user_util = (cpu.user_time / wall_time) * 100.0;
  double kernel_util = (cpu.kernel_time / wall_time) * 100.0;

  std::cout << "-----------------------------------------------" << std::endl;
  printf("Total Computation Time:    %.3f seconds\n", computation_time);
  printf("Start-to-End Wall Time:    %.3f seconds\n", wall_time);
  std::cout << std::endl;
  printf("CPU Utilization:           %.2f %%  +  %.2f %% kernel overhead\n",
         user_util, kernel_util);
  printf("Multi-core Efficiency:     %.2f %%  +  %.2f %% kernel overhead\n",
         user_util / threads, kernel_util / threads);
  std::cout << "-----------------------------------------------" << std::endl;

  // Generate the professional validation report
  write_validation_report(digits, computation_time, wall_time, user_util, kernel_util, event_history, result_str, actual_len);

  delete[] result_str;
  mpz_clears(pi_z, num, sqrt_val, d10, NULL);
  return 0;
}
