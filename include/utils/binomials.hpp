#pragma once
#include <NTL/RR.h>
#include <NTL/ZZ.h>
#include <NTL/matrix.h>
#include <cstdint>
#include <cmath>

/* binomials are precomputed up to MAX_N-choose-MAX_T */
#define MAX_N 2000
#define MAX_T 300

#define LOW_K_MAX_N 10000
#define LOW_K_MAX_T 10

extern NTL::RR pi;
extern NTL::RR nat_log_2;

extern NTL::Mat<NTL::ZZ> binomial_table;
extern NTL::Mat<NTL::ZZ> low_k_binomial_table;
extern bool is_data_initialized;

void InitBinomials();

NTL::RR lnFactorial(NTL::RR n);
NTL::RR lnBinom(NTL::RR n, NTL::RR k);
NTL::ZZ binomial_wrapper(long n, long k);
NTL::RR log2_RR(NTL::RR v);
