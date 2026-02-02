#pragma once
#include <algorithm>
#include <stdint.h>
#include <string>
#include <vector>


namespace pi {

class BigInt {
public:
  std::vector<uint64_t> data;
  bool negative = false;

  BigInt() {}
  BigInt(int64_t val) {
    if (val < 0) {
      negative = true;
      val = -val;
    }
    if (val)
      data.push_back((uint64_t)val);
  }

  bool is_zero() const { return data.empty(); }

  void add(const BigInt &other);
  void sub(const BigInt &other);
  void mul(const BigInt &other);
  void mul_small(uint64_t val);

  static void binary_split(int64_t a, int64_t b, BigInt &P, BigInt &Q,
                           BigInt &T);

  void shift_left(size_t limbs);
  void shift_right(size_t limbs);

  static BigInt sqrt(const BigInt &s, int64_t precision_limbs);
  static BigInt div(const BigInt &a, const BigInt &b, int64_t precision_limbs);

  std::string to_string() const;
  std::string divide_to_string(const BigInt &divisor, int precision) const;
};

} // namespace pi
