#pragma once
#include <stdint.h>
#include <vector>


namespace pi {

class NTT {
public:
  struct Prime {
    uint64_t p;
    uint64_t g;
  };

  static constexpr Prime P1 = {998244353, 3};
  static constexpr Prime P2 = {1004535809, 3};
  static constexpr Prime P3 = {469762049, 3};

  static uint64_t power(uint64_t a, uint64_t b, uint64_t p);
  static uint64_t modInverse(uint64_t n, uint64_t p);
  static void ntt(std::vector<uint64_t> &a, bool invert, Prime prime);

  static std::vector<uint64_t> multiply(const std::vector<uint64_t> &a,
                                        const std::vector<uint64_t> &b);

private:
  static std::vector<uint64_t> ntt_single_prime(const std::vector<uint64_t> &a,
                                                const std::vector<uint64_t> &b,
                                                Prime prime);
};

} // namespace pi
