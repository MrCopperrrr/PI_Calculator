#pragma once
#include <cstdint>
#include <gmp.h>
#include <string>
#include <vector>


namespace pi {

class BaseConverter {
public:
  static void parallel_to_str(mpz_t n, int64_t total_digits, char *out_buf);

private:
  static void recursive_split(mpz_t n, int64_t digits, char *out,
                              const std::vector<mpz_t *> &powers, int level);
};

} // namespace pi
