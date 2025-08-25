#pragma once
#include "pqt.h"

// This is now a serial binary split function, not parallel.
PQT binary_split(long long a, long long b, const mpz_class& C3_over_24);