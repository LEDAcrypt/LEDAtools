#pragma once
#include "binomials.hpp"
#include "proper_primes.hpp"
#include <NTL/RR.h>
#include <utility>

// choice of the approximation praxis for the estimated fraction of an error
// to appear in the next iteration of a bit-flipping decoder
#define ROUNDING_PRAXIS round

NTL::RR compute_p_cc(const uint64_t d_c, const uint64_t n, const uint64_t t);
NTL::RR compute_p_ci(const uint64_t d_c, const uint64_t n, const uint64_t t);
NTL::RR compute_p_ic(const uint64_t d_c, const uint64_t n, const uint64_t t);
NTL::RR compute_p_ii(const uint64_t d_c, const uint64_t n, const uint64_t t);
NTL::RR ComputePrBitCorrection(const NTL::RR p_ic, const uint64_t d_v,
                               const uint64_t threshold);
NTL::RR ComputePrBitFaultInduction(const NTL::RR p_ci, const uint64_t d_v,
                                   const uint64_t threshold);
NTL::RR ComputePrBitCorrectionMulti(const NTL::RR p_ic, const uint64_t d_v,
                                    const uint64_t t, const uint64_t threshold,
                                    const uint64_t toCorrect);
NTL::RR ComputePrBitInduceMulti(const NTL::RR p_ci, const uint64_t d_v,
                                const uint64_t t, const uint64_t n,
                                const uint64_t threshold,
                                const uint64_t toInduce);
uint64_t FindNextNumErrors(const uint64_t n_0, const uint64_t p,
                           const uint64_t d_v, const uint64_t t);
std::pair<NTL::RR, uint64_t> Find1IterDFR(const uint64_t n_0, const uint64_t p,
                                          const uint64_t d_v, const uint64_t t);
std::pair<NTL::RR, uint64_t>
Find1IterTLeftoverPr(const uint64_t n_0, const uint64_t p, const uint64_t d_v,
                     const uint64_t t, const uint64_t t_leftover);
uint64_t Findpth(const uint64_t n_0, const uint64_t d_v_prime,
                 const uint64_t t);
