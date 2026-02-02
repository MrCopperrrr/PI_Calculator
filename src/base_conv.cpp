#include "base_conv.hpp"
#include <cmath>
#include <cstring>
#include <iostream>
#include <omp.h>


namespace pi {

struct MpzWrapper {
  mpz_t val;
  MpzWrapper() { mpz_init(val); }
  ~MpzWrapper() { mpz_clear(val); }
};

void BaseConverter::recursive_split(mpz_t n, int64_t digits, char *out,
                                    const std::vector<mpz_t *> &powers,
                                    int level) {
  if (digits <= 16384) {
    char *s = mpz_get_str(NULL, 10, n);
    size_t len = strlen(s);

    if ((int64_t)len < digits) {
      memset(out, '0', digits - len);
    }
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

#pragma omp task shared(out, powers) firstprivate(high, digits, half, level)
  {
    mpz_t h;
    mpz_init_set(h, high);
    recursive_split(h, digits - half, out, powers, level + 1);
    mpz_clear(h);
  }

#pragma omp task shared(out, powers) firstprivate(low, digits, half, level)
  {
    mpz_t l;
    mpz_init_set(l, low);
    recursive_split(l, half, out + (digits - half), powers, level + 1);
    mpz_clear(l);
  }

#pragma omp taskwait
  mpz_clear(high);
  mpz_clear(low);
}

void BaseConverter::parallel_to_str(mpz_t n, int64_t total_digits,
                                    char *out_buf) {
  std::vector<int64_t> d_steps;
  int64_t cur_d = total_digits;
  while (cur_d > 16384) {
    int64_t half = cur_d / 2;
    d_steps.push_back(half);
    cur_d = cur_d - half;
  }

  std::vector<MpzWrapper *> wrappers(d_steps.size());
  std::vector<mpz_t *> powers(d_steps.size());

#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < (int)d_steps.size(); ++i) {
    wrappers[i] = new MpzWrapper();
    mpz_ui_pow_ui(wrappers[i]->val, 10, d_steps[i]);
    powers[i] = &(wrappers[i]->val);
  }

#pragma omp parallel
  {
#pragma omp single
    recursive_split(n, total_digits, out_buf, powers, 0);
  }

  out_buf[total_digits] = '\0';
  for (auto w : wrappers)
    delete w;
}

} // namespace pi
