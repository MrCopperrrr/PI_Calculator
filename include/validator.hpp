#ifndef VALIDATOR_HPP
#define VALIDATOR_HPP

#include <cstdint>
#include <string>
#include <vector>
#include <utility>

namespace pi {

struct ValidationResult {
  bool is_known_benchmark;
  bool spot_check_passed;
  std::string expected_last_digits;
  std::string actual_last_digits;
  std::vector<int64_t> digit_counts; // Counts of '0' through '9'
  uint64_t dec_hash;                 // Modulo 2^61 - 1 hash
  double chi_square;                 // Statistical uniformity check
};

class PiValidator {
public:
  // Validate computed Pi decimal string against known benchmarks and statistics
  static ValidationResult validate(const char *pi_str, int64_t total_digits);

  // Generate y-cruncher style validation file and console summary
  static void write_validation_file(
      const char *filename, const char *pi_str, int64_t total_digits,
      const ValidationResult &val_res, double comp_time, double wall_time,
      double user_time, double kernel_time, int threads,
      const std::vector<std::pair<double, std::string>> &event_log);

private:
  static uint64_t calculate_hash(const char *str, int64_t len);
};

} // namespace pi

#endif // VALIDATOR_HPP
