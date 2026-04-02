#pragma once
#include <cstdint>
#include <gmp.h>
#include <iostream>
#include <string>

namespace pi {

class BigInt {
public:
  mpz_t value;

  BigInt() { mpz_init(value); }
  BigInt(int64_t val) { mpz_init_set_si(value, val); }

  BigInt(const BigInt &other) { mpz_init_set(value, other.value); }
  BigInt(BigInt &&other) noexcept {
    value[0] = other.value[0];
    mpz_init(other.value);
  }
  BigInt &operator=(const BigInt &other) {
    if (this != &other)
      mpz_set(value, other.value);
    return *this;
  }
  BigInt &operator=(BigInt &&other) noexcept {
    if (this != &other) {
      mpz_clear(value);
      value[0] = other.value[0];
      mpz_init(other.value);
    }
    return *this;
  }

  ~BigInt() { mpz_clear(value); }

  bool is_zero() const { return mpz_sgn(value) == 0; }

  void add(const BigInt &other) { mpz_add(value, value, other.value); }
  void sub(const BigInt &other) { mpz_sub(value, value, other.value); }
  void mul(const BigInt &other) { mpz_mul(value, value, other.value); }
  void mul_small(uint64_t val) { mpz_mul_ui(value, value, val); }

  static void binary_split(int64_t a, int64_t b, BigInt &P, BigInt &Q,
                           BigInt &T);

  static void parallel_pow_ui(mpz_t rop, uint64_t base, uint64_t exp);
  static void parallel_sqrt(mpz_t rop, const mpz_t n);
  static void parallel_div(mpz_t q, const mpz_t num, const mpz_t den);

  void shift_left(size_t limbs) { mpz_mul_2exp(value, value, limbs * 64); }

  std::string to_string() const {
    char *buf = mpz_get_str(NULL, 10, value);
    std::string s(buf);
    void (*freefunc)(void *, size_t);
    mp_get_memory_functions(NULL, NULL, &freefunc);
    freefunc(buf, s.size() + 1);
    return s;
  }
};

} // namespace pi
