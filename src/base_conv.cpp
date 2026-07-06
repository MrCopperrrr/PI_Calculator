#include "base_conv.hpp"
#include "bigint.hpp"
#include "ntt.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <omp.h>
#include <vector>

namespace pi {

struct PowerPair {
  int64_t digits;
  mpz_t val;
  PowerPair() { mpz_init(val); }
  PowerPair(const PowerPair &other) {
    digits = other.digits;
    mpz_init_set(val, other.val);
  }
  PowerPair(PowerPair &&other) noexcept {
    digits = other.digits;
    val[0] = other.val[0];
    mpz_init(other.val);
  }
  PowerPair &operator=(const PowerPair &other) {
    if (this != &other) {
      digits = other.digits;
      mpz_set(val, other.val);
    }
    return *this;
  }
  PowerPair &operator=(PowerPair &&other) noexcept {
    if (this != &other) {
      digits = other.digits;
      mpz_clear(val);
      val[0] = other.val[0];
      mpz_init(other.val);
    }
    return *this;
  }
  ~PowerPair() { mpz_clear(val); }
};

void BaseConverter::collect_powers(int64_t digits, std::vector<int64_t> &needed) {
  if (digits <= 16384)
    return;
  int64_t half = digits / 2;
  needed.push_back(half);
  collect_powers(digits - half, needed);
  collect_powers(half, needed);
}

static const mpz_t *get_power(const std::vector<PowerPair> &powers, int64_t half) {
  int l = 0, r = (int)powers.size() - 1;
  while (l <= r) {
    int m = l + (r - l) / 2;
    if (powers[m].digits == half)
      return &powers[m].val;
    if (powers[m].digits < half)
      l = m + 1;
    else
      r = m - 1;
  }
  return nullptr;
}

void BaseConverter::recursive_split(mpz_t n, int64_t digits, char *out,
                                    const std::vector<mpz_t *> &powers,
                                    int level) {
  // Use a much higher threshold for tasking to avoid memory bloat
  // 1 million digits is a good balance between parallelism and memory safety
  if (digits < 1000000) {
    if (digits <= 16384) {
      char *s = mpz_get_str(NULL, 10, n);
      size_t len = strlen(s);
      if ((int64_t)len < digits) memset(out, '0', digits - len);
      memcpy(out + (digits > (int64_t)len ? digits - len : 0), s, len);
      void (*freefunc)(void *, size_t);
      mp_get_memory_functions(NULL, NULL, &freefunc);
      freefunc(s, len + 1);
      return;
    }

    int64_t half = digits / 2;
    mpz_t high, low;
    mpz_init(high);
    mpz_init(low);
    mpz_tdiv_qr(high, low, n, *powers[level]);
    recursive_split(high, digits - half, out, powers, level + 1);
    recursive_split(low, half, out + (digits - half), powers, level + 1);
    mpz_clear(high);
    mpz_clear(low);
    return;
  }

  int64_t half = digits / 2;
  const mpz_t *power = get_power(powers, half);
  mpz_t high, low;
  mpz_init(high);
  mpz_init(low);
  mpz_tdiv_qr(high, low, n, *powers[level]);

#pragma omp task shared(out, powers, high) firstprivate(digits, half)
  {
    mpz_t h;
    mpz_init_set(h, high);
    recursive_split(h, digits - half, out, powers, level + 1);
    mpz_clear(h);
    mpz_clear(high); // Clear the firstprivate copy
  }

#pragma omp task shared(out, powers, low) firstprivate(digits, half)
  {
    mpz_t l;
    mpz_init_set(l, low);
    recursive_split(l, half, out + (digits - half), powers, level + 1);
    mpz_clear(l);
    mpz_clear(low); // Clear the firstprivate copy
  }

#pragma omp taskwait
}

void BaseConverter::parallel_to_str(mpz_t n, int64_t total_digits,
                                    char *out_buf) {
  std::vector<int64_t> needed;
  collect_powers(total_digits, needed);
  std::sort(needed.begin(), needed.end());
  needed.erase(std::unique(needed.begin(), needed.end()), needed.end());

  std::vector<PowerPair> powers(needed.size());
  for (size_t i = 0; i < needed.size(); ++i) {
    powers[i].digits = needed[i];
    bool computed = false;
    for (int j = (int)i - 1; j >= 0; --j) {
      if (powers[j].digits * 2 == powers[i].digits) {
        NTTMultiplier::multiply(powers[i].val, powers[j].val, powers[j].val);
        computed = true;
        break;
      } else if (powers[j].digits < powers[i].digits) {
        int64_t rem = powers[i].digits - powers[j].digits;
        for (int k = j; k >= 0; --k) {
          if (powers[k].digits == rem) {
            NTTMultiplier::multiply(powers[i].val, powers[j].val,
                                    powers[k].val);
            computed = true;
            break;
          }
        }
        if (computed)
          break;
      }
    }
    if (!computed) {
      BigInt::parallel_pow_ui(powers[i].val, 10, powers[i].digits);
    }
  }

#pragma omp parallel
  {
#pragma omp single
    recursive_split(n, total_digits, out_buf, powers);
  }

  out_buf[total_digits] = '\0';
}

} // namespace pi
