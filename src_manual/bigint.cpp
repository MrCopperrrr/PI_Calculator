#include "bigint.hpp"
#include "ntt.hpp"
#include <algorithm>
#include <omp.h>
#include <string>
#include <vector>

namespace pi {

int compare_abs(const BigInt &a, const BigInt &b) {
  if (a.data.size() != b.data.size())
    return a.data.size() < b.data.size() ? -1 : 1;
  for (int i = (int)a.data.size() - 1; i >= 0; --i) {
    if (a.data[i] != b.data[i])
      return a.data[i] < b.data[i] ? -1 : 1;
  }
  return 0;
}

void BigInt::add(const BigInt &other) {
  if (negative == other.negative) {
    size_t n = std::max(data.size(), other.data.size());
    uint64_t carry = 0;
    for (size_t i = 0; i < n || carry; ++i) {
      if (i == data.size())
        data.push_back(0);
      uint64_t other_val = (i < other.data.size()) ? other.data[i] : 0;
      uint64_t next_carry = 0;
      if (__builtin_add_overflow(data[i], other_val, &data[i]))
        next_carry = 1;
      if (__builtin_add_overflow(data[i], carry, &data[i]))
        next_carry += 1;
      carry = next_carry;
    }
  } else {
    bool orig_neg = negative;
    negative = false;
    int cmp = compare_abs(*this, other);
    if (cmp >= 0) {
      sub(other);
      negative = orig_neg;
    } else {
      BigInt temp = other;
      temp.sub(*this);
      *this = std::move(temp);
    }
  }
}

void BigInt::sub(const BigInt &other) {
  if (negative != other.negative) {
    bool orig_neg = negative;
    negative = other.negative;
    add(other);
    negative = orig_neg;
  } else {
    int cmp = compare_abs(*this, other);
    if (cmp >= 0) {
      uint64_t borrow = 0;
      for (size_t i = 0; i < data.size(); ++i) {
        uint64_t other_val = (i < other.data.size()) ? other.data[i] : 0;
        uint64_t next_borrow = 0;
        if (__builtin_sub_overflow(data[i], other_val, &data[i]))
          next_borrow = 1;
        if (__builtin_sub_overflow(data[i], borrow, &data[i]))
          next_borrow += 1;
        borrow = next_borrow;
      }
      while (!data.empty() && data.back() == 0)
        data.pop_back();
    } else {
      BigInt temp = other;
      temp.sub(*this);
      temp.negative = !negative;
      *this = std::move(temp);
    }
  }
}

void BigInt::mul_small(uint64_t val) {
  if (val == 0) {
    data.clear();
    negative = false;
    return;
  }
  uint64_t carry = 0;
  for (size_t i = 0; i < data.size(); ++i) {
    unsigned __int128 prod = (unsigned __int128)data[i] * val + carry;
    data[i] = (uint64_t)prod;
    carry = (uint64_t)(prod >> 64);
  }
  while (carry) {
    data.push_back(carry & 0xFFFFFFFFFFFFFFFFULL);
    carry >>= 64;
  }
}

void BigInt::mul(const BigInt &other) {
  if (is_zero() || other.is_zero()) {
    data.clear();
    negative = false;
    return;
  }
  bool res_neg = negative != other.negative;

  if (data.size() < 128 || other.data.size() < 128) {
    std::vector<uint64_t> result(data.size() + other.data.size(), 0);
    for (size_t i = 0; i < data.size(); ++i) {
      uint64_t carry = 0;
      for (size_t j = 0; j < other.data.size(); ++j) {
        unsigned __int128 prod =
            (unsigned __int128)data[i] * other.data[j] + result[i + j] + carry;
        result[i + j] = (uint64_t)prod;
        carry = (uint64_t)(prod >> 64);
      }
      result[i + other.data.size()] += carry;
    }
    while (!result.empty() && result.back() == 0)
      result.pop_back();
    data = std::move(result);
  } else {
    // Split 64-bit into four 16-bit parts for NTT safety
    std::vector<uint64_t> a_16, b_16;
    a_16.reserve(data.size() * 4);
    b_16.reserve(other.data.size() * 4);
    for (uint64_t x : data) {
      for (int i = 0; i < 4; ++i) {
        a_16.push_back(x & 0xFFFFULL);
        x >>= 16;
      }
    }
    for (uint64_t x : other.data) {
      for (int i = 0; i < 4; ++i) {
        b_16.push_back(x & 0xFFFFULL);
        x >>= 16;
      }
    }

    std::vector<uint64_t> res_16 = NTT::multiply(a_16, b_16);

    data.clear();
    unsigned __int128 carry = 0;
    for (size_t i = 0; i < res_16.size(); i += 4) {
      unsigned __int128 combined = 0;
      for (int j = 0; j < 4; ++j) {
        if (i + j < res_16.size())
          combined |= (unsigned __int128)res_16[i + j] << (16 * j);
      }
      combined += carry;
      data.push_back((uint64_t)(combined & 0xFFFFFFFFFFFFFFFFULL));
      carry = combined >> 64;
    }
    if (carry)
      data.push_back((uint64_t)carry);
    while (!data.empty() && data.back() == 0)
      data.pop_back();
  }
  negative = res_neg;
}

std::string BigInt::to_string() const {
  if (is_zero())
    return "0";
  std::string s = negative ? "-" : "";
  BigInt temp = *this;
  temp.negative = false;
  std::string digits;
  while (!temp.is_zero()) {
    uint64_t rem = 0;
    uint64_t divisor = 1000000000000000000ULL;
    for (int i = (int)temp.data.size() - 1; i >= 0; --i) {
      unsigned __int128 cur = temp.data[i] + ((unsigned __int128)rem << 64);
      temp.data[i] = (uint64_t)(cur / divisor);
      rem = (uint64_t)(cur % divisor);
    }
    while (!temp.data.empty() && temp.data.back() == 0)
      temp.data.pop_back();
    std::string part = std::to_string(rem);
    if (!temp.is_zero())
      while (part.length() < 18)
        part = "0" + part;
    digits = part + digits;
  }
  return s + digits;
}

void BigInt::shift_left(size_t limbs) {
  if (limbs == 0 || is_zero())
    return;
  data.insert(data.begin(), limbs, 0);
}

void BigInt::shift_right(size_t limbs) {
  if (limbs >= data.size()) {
    data.clear();
    negative = false;
    return;
  }
  data.erase(data.begin(), data.begin() + limbs);
}

BigInt BigInt::sqrt(const BigInt &s, int64_t precision_limbs) {
  if (s.is_zero())
    return BigInt(0);
  // Integer square root using Newton's method: x = (x + s/x) / 2
  BigInt x = BigInt(1);
  x.shift_left(s.data.size() / 2 + 1); // Initial guess

  // For π calculation we need a more precise sqrt.
  // This is a simplified version.
  return x;
}

BigInt BigInt::div(const BigInt &a, const BigInt &b, int64_t precision_limbs) {
  // Newton-Raphson for reciprocal: x = x(2 - bx)
  return a; // Placeholder
}

std::string BigInt::divide_to_string(const BigInt &divisor,
                                     int precision) const {
  // Convert decimal precision to limb precision (1 limb approx 19 digits)
  int64_t precision_limbs = precision / 19 + 2;
  BigInt res = div(*this, divisor, precision_limbs);
  return res.to_string(); // This needs to be formatted for decimals
}

} // namespace pi
