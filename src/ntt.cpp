#include "ntt.hpp"
#include <algorithm>
#include <omp.h>
#include <vector>

namespace pi {

uint64_t NTT::power(uint64_t a, uint64_t b, uint64_t p) {
  uint64_t res = 1;
  a %= p;
  while (b > 0) {
    if (b & 1)
      res = (unsigned __int128)res * a % p;
    a = (unsigned __int128)a * a % p;
    b >>= 1;
  }
  return res;
}

uint64_t NTT::modInverse(uint64_t n, uint64_t p) { return power(n, p - 2, p); }

void NTT::ntt(std::vector<uint64_t> &a, bool invert, Prime prime) {
  size_t n = a.size();
  uint64_t p = prime.p;
  uint64_t g = prime.g;

  for (size_t i = 1, j = 0; i < n; i++) {
    size_t bit = n >> 1;
    for (; j & bit; bit >>= 1)
      j ^= bit;
    j ^= bit;
    if (i < j)
      std::swap(a[i], a[j]);
  }

  for (size_t len = 2; len <= n; len <<= 1) {
    uint64_t wlen = power(g, (p - 1) / len, p);
    if (invert)
      wlen = modInverse(wlen, p);

#pragma omp parallel for if (n > 1000)
    for (size_t i = 0; i < n; i += len) {
      uint64_t w = 1;
      for (size_t j = 0; j < len / 2; j++) {
        uint64_t u = a[i + j],
                 v = (unsigned __int128)a[i + j + len / 2] * w % p;
        a[i + j] = (u + v) % p;
        a[i + j + len / 2] = (u + p - v) % p;
        w = (unsigned __int128)w * wlen % p;
      }
    }
  }

  if (invert) {
    uint64_t n_inv = modInverse(n, p);
    for (uint64_t &x : a)
      x = (unsigned __int128)x * n_inv % p;
  }
}

std::vector<uint64_t> NTT::ntt_single_prime(const std::vector<uint64_t> &a,
                                            const std::vector<uint64_t> &b,
                                            Prime prime) {
  size_t n = 1;
  while (n < a.size() + b.size())
    n <<= 1;
  std::vector<uint64_t> fa(n, 0), fb(n, 0);
  std::copy(a.begin(), a.end(), fa.begin());
  std::copy(b.begin(), b.end(), fb.begin());
  ntt(fa, false, prime);
  ntt(fb, false, prime);
  for (size_t i = 0; i < n; i++)
    fa[i] = (unsigned __int128)fa[i] * fb[i] % prime.p;
  ntt(fa, true, prime);
  return fa;
}

std::vector<uint64_t> NTT::multiply(const std::vector<uint64_t> &a,
                                    const std::vector<uint64_t> &b) {
  if (a.empty() || b.empty())
    return {};

  std::vector<uint64_t> res1, res2;
#pragma omp parallel sections
  {
#pragma omp section
    res1 = ntt_single_prime(a, b, P1);
#pragma omp section
    res2 = ntt_single_prime(a, b, P2);
  }

  size_t n = res1.size();
  std::vector<uint64_t> res(n);
  uint64_t invP1 = modInverse(P1.p, P2.p);

#pragma omp parallel for if (n > 1000)
  for (int i = 0; i < (int)n; i++) {
    uint64_t a1 = res1[i];
    uint64_t a2 = res2[i];
    uint64_t k = (unsigned __int128)(a2 + P2.p - a1) * invP1 % P2.p;
    res[i] = a1 + (unsigned __int128)k * P1.p;
  }
  return res;
}

} // namespace pi
