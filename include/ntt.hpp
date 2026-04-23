#pragma once
#include <cstdint>
#include <gmp.h>
#include <omp.h>
#include <vector>


namespace pi {

class NTTMultiplier {
public:
  // Primes for Triple-Prime NTT to ensure accuracy for huge numbers
  static constexpr uint64_t MODS[] = {998244353, 1004535809, 469762049};
  static constexpr uint64_t G = 3;

  static void ntt(std::vector<uint64_t> &a, bool invert, uint64_t mod);

  // Multiplies two mpz_t using parallel NTT
  static void multiply(mpz_t rop, const mpz_t op1, const mpz_t op2);
  
  static bool use_hybrid; // Flag to toggle between Step 1 and Step 2 strategies

private:
  static uint64_t power(uint64_t base, uint64_t exp, uint64_t mod);
  static uint64_t modInverse(uint64_t n, uint64_t mod);

  // Helpers to convert mpz_t to/from NTT buffers
  static std::vector<uint64_t> mpz_to_vec(const mpz_t n, size_t &limbs);
  static void vec_to_mpz(mpz_t rop, const std::vector<uint64_t> &vec);
};

} // namespace pi
