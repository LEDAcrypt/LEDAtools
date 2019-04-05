#pragma once
#include <NTL/RR.h>
#include <utility>
#include "proper_primes.hpp"

// choice of the approximation praxis for the estimated fraction of an error
// to appear in the next iteration of a bit-flipping decoder
#define ROUNDING_PRAXIS round

/* Probability that a variable node is correct, and a parity equation involving
 * it is satisfied */
NTL::RR compute_p_cc(const uint64_t d_c, 
                   const uint64_t n, 
                   const uint64_t t){
    NTL::RR result = NTL::RR(0);
    uint64_t bound = (d_c - 1) < t ? d_c - 1 : t;

    /* the number of errors falling in the PC equation should be at least 
     * the amount which cannot be placed in a non checked place          */
    uint64_t LowerTHitBound = (n-d_c) < t ? t-(n-d_c) : 0;
    /* and it should be even, since the PC equation must be satisfied */
    LowerTHitBound = LowerTHitBound % 2 ? LowerTHitBound + 1 : LowerTHitBound;

    for(uint64_t j = LowerTHitBound; j <= bound; j = j+2 ){
        result += to_RR( binomial_wrapper(d_c-1,j) * binomial_wrapper(n-d_c,t-j) ) / 
                         to_RR( binomial_wrapper(n-1,t) );
    }
    return result;
}

/* Probability that a variable node is correct, and a parity equation involving
 * it is *not* satisfied */
NTL::RR compute_p_ci(const uint64_t d_c, 
                   const uint64_t n,
                   const uint64_t t){
    NTL::RR result = NTL::RR(0);
    uint64_t bound = (d_c - 1) < t ? d_c - 1 : t;

    /* the number of errors falling in the PC equation should be at least 
     * the amount which cannot be placed in a non checked place          */
    uint64_t LowerTHitBound = (n-d_c) < t ? t-(n-d_c) : 1;
    /* and it should be odd, since the PC equation must be non satisfied */
    LowerTHitBound = LowerTHitBound % 2 ? LowerTHitBound : LowerTHitBound + 1;

    for(uint64_t j = LowerTHitBound; j <= bound; j = j+2 ){
        result += to_RR( binomial_wrapper(d_c-1,j) * binomial_wrapper(n-d_c,t-j) )
                 / to_RR( binomial_wrapper(n-1,t) );
    }
    return result;
}

/* Probability that a variable node is *not* correct, and a parity equation involving
 * it is *not* satisfied */
NTL::RR compute_p_ic(const uint64_t d_c, 
                   const uint64_t n,
                   const uint64_t t){
    NTL::RR result = NTL::RR(0);
    uint64_t UpperTBound = (d_c - 1) < t - 1 ? d_c - 1 : t - 1;

    /* the number of errors falling in the PC equation should be at least 
     * the amount which cannot be placed in a non checked place          */
    uint64_t LowerTHitBound = (n-d_c-1) < (t-1) ? (t-1)-(n-d_c-1) : 0;
    /* and it should be even, since the PC equation must be unsatisfied (when 
     * accounting for the one we are considering as already placed*/
    LowerTHitBound = LowerTHitBound % 2 ? LowerTHitBound + 1 : LowerTHitBound;

    for(uint64_t j = LowerTHitBound; j <= UpperTBound; j = j+2 ){
        result += NTL::to_RR( binomial_wrapper(d_c-1,j) * binomial_wrapper(n-d_c,t-j-1) ) 
                 / to_RR( binomial_wrapper(n-1,t-1) );
    }
    return result;
}

/* Probability that a variable node is *not* correct, and a parity equation involving
 * it is satisfied */
NTL::RR compute_p_ii(const uint64_t d_c, 
                   const uint64_t n,
                   const uint64_t t){

    NTL::RR result = NTL::RR(0);
    uint64_t bound = (d_c - 1) < t - 1 ? d_c - 1 : t - 1;
    
    /* the number of errors falling in the PC equation should be at least 
     * the amount which cannot be placed in a non checked place          */
    uint64_t LowerTHitBound = (n-d_c) < (t-1) ? (t-1)-(n-d_c) : 1;
    /* and it should be odd, since the PC equation must be satisfied (when 
     * accounting for the one we are considering as already placed)*/
    LowerTHitBound = LowerTHitBound % 2 ? LowerTHitBound : LowerTHitBound +1;
    for(uint64_t j = LowerTHitBound; j <= bound; j = j+2 ){
        result += NTL::to_RR( binomial_wrapper(d_c-1,j) * binomial_wrapper(n-d_c,t-j-1) ) 
                 / to_RR( binomial_wrapper(n-1,t-1) );
    }
    return result;
}

/* note p_cc + p_ci = 1 */
/* note p_ic + p_ii = 1 */

/* Probability that a given erroneous variable is deemed as such, and is thus
 * corrected, given a threshold for the amount of unsatisfied parity check
 * equations. Called P_ic in most texts */
NTL::RR ComputePrBitCorrection( const NTL::RR p_ic, 
                                const uint64_t d_v,
                                const uint64_t t,
                                const uint64_t threshold ){
// 		Pic=0; /* p_correct */
// 		for (j=b,dv,
// 			term=binomial(dv,j)*(p_ic^j)*(1-p_ic)^(dv-j);
// 			Pic=Pic+term;
// 		);
  NTL::RR result = NTL::RR(0), success, failure;
  for (uint64_t j = threshold; j <= d_v; j++){
     NTL::pow(success, p_ic, NTL::to_RR(j));
     NTL::pow(failure, NTL::RR(1)-p_ic, NTL::to_RR(d_v-j));
     result += NTL::to_RR(binomial_wrapper(d_v,j)) * success * failure;
  }
  return result;
}

/* Probability that a given correct variable is not deemed as such, and is thus
 * fault-induced, given a threshold for the amount of unsatisfied parity check
 * equations. Called P_ci in most texts, p_induce in official comment */
NTL::RR ComputePrBitFaultInduction( const NTL::RR p_ci,
                                    const uint64_t d_v,
                                    const uint64_t t, /* unused */
                                    const uint64_t threshold ){

  NTL::RR result= NTL::RR(0), success, failure;
  for (uint64_t j = threshold; j <= d_v; j++){
     NTL::pow(success, p_ci, NTL::to_RR(j));
     NTL::pow(failure, NTL::RR(1)-p_ci, NTL::to_RR(d_v-j));
     result += NTL::to_RR(binomial_wrapper(d_v,j)) * success * failure;
  }
  return result;
}

/* computes the probability that toCorrect bits are corrected
 * known as P{N_ic = toCorrect}  */
NTL::RR ComputePrBitCorrectionMulti( const NTL::RR p_ic, 
                                const uint64_t d_v,
                                const uint64_t t,
                                const uint64_t threshold,
                                const uint64_t toCorrect){
   NTL::RR ProbCorrectOne = ComputePrBitCorrection(p_ic,d_v,t,threshold);
   return NTL::to_RR(binomial_wrapper(t,toCorrect)) * 
          NTL::pow(ProbCorrectOne,NTL::RR(toCorrect)) *
          NTL::pow(1-ProbCorrectOne,NTL::RR(t-toCorrect));
}

/* computes the probability that toInduce faults are induced 
 * known as P{N_ci = toInduce} or Pr{f_wrong = to_induce} */
NTL::RR ComputePrBitInduceMulti(const NTL::RR p_ci, 
                                const uint64_t d_v,
                                const uint64_t t,
                                const uint64_t n,
                                const uint64_t threshold,
                                const uint64_t toInduce){
//    if(toInduce <= 1 ){
//        return NTL::RR(0);
//    }    
   NTL::RR ProbInduceOne = ComputePrBitFaultInduction(p_ci,d_v,t,threshold);
   return NTL::to_RR(binomial_wrapper(n-t,toInduce)) * 
          NTL::pow(ProbInduceOne,NTL::RR(toInduce)) *
          NTL::pow(1-ProbInduceOne,NTL::RR(n-t-toInduce));                                    
}

uint64_t FindNextNumErrors(const uint64_t n_0,
                           const uint64_t p,
                           const uint64_t d_v,
                           const uint64_t t){
    NTL::RR p_ci, p_ic;
     p_ci = compute_p_ci(n_0*d_v,n_0*p,t);
     p_ic = compute_p_ic(n_0*d_v,n_0*p,t);
    uint64_t t_next=t;
//      uint64_t best_threshold = (d_v - 1)/2;
    for(uint64_t i = (d_v - 1)/2; i <= d_v - 1; i++){
       NTL::RR t_approx=  t -
                          t * ComputePrBitCorrection(p_ic, d_v, t, i) +
                          (n_0*p - t) * ComputePrBitFaultInduction(p_ci, d_v, t, i);
       unsigned long int t_curr = NTL::conv<unsigned long int>(NTL::ROUNDING_PRAXIS(t_approx)) ;
       /*Note : we increase the threshold only if it improves strictly on the 
        * predicted error correction. */
       if (t_curr < t_next){
          t_next = t_curr;
//           best_threshold = i;
       }
    }
    /* considering that any code will correct a single bit error, if 
     * t_next == 1, we save a computation iteration and shortcut to t_next == 0*/
    if (t_next == 1) {
        t_next = 0;
    }
    return t_next;
}

/* computes the exact 1-iteration DFR and the best threshold on the number of
 * upcs to achieve it */
std::pair<NTL::RR,uint64_t> Find1IterDFR(const uint64_t n_0,
                                         const uint64_t p,
                                         const uint64_t d_v,
                                         const uint64_t t){
    NTL::RR p_ci, p_ic, P_correct, P_induce;
    NTL::RR DFR, best_DFR = NTL::RR(1);
    p_ci = compute_p_ci(n_0*d_v,n_0*p,t);
    p_ic = compute_p_ic(n_0*d_v,n_0*p,t);
    uint64_t best_threshold = (d_v - 1)/2;
    for(uint64_t b = best_threshold; b <= d_v - 1; b++){
       DFR = NTL::RR(1) - ComputePrBitCorrectionMulti(p_ic, d_v, t, b, t) * ComputePrBitInduceMulti(p_ci,d_v,t,n_0*p,b,0);
       /*Note : we increase the threshold only if it improves strictly on the 
        * predicted error correction. */
       if (DFR < best_DFR){
          best_DFR = DFR;
          best_threshold = b;
       }
    }
//     std::cout << best_threshold << std::endl;
    return std::make_pair(best_DFR,best_threshold);
}


/* computes the exact 1-iteration probability of leaving at most t_leftover
 * uncorrected errors out of t. */
std::pair<NTL::RR,uint64_t> Find1IterTLeftoverPr(const uint64_t n_0,
                                         const uint64_t p,
                                         const uint64_t d_v,
                                         const uint64_t t,
                                         const uint64_t t_leftover){
    NTL::RR p_ci, p_ic;
    NTL::RR DFR, best_DFR = NTL::RR(1);
    p_ci = compute_p_ci(n_0*d_v,n_0*p,t);
    p_ic = compute_p_ic(n_0*d_v,n_0*p,t);
    int n= p*n_0;
    uint64_t best_threshold = (d_v + 1)/2;
    
    for(uint64_t b = best_threshold; b <= d_v ; b++){
       DFR = NTL::RR(0);
       NTL::RR P_correct = ComputePrBitCorrection(p_ic, d_v, t,b);
       NTL::RR P_induce = ComputePrBitFaultInduction(p_ci,d_v, t/* unused */,b);
       for(int tau = 0 ; tau <= t_leftover; tau++){
         for(int n_to_induce = 0 ; n_to_induce <= t_leftover; n_to_induce++) {
             NTL::RR prob_induce_n = NTL::to_RR(binomial_wrapper(n-t,n_to_induce)) *
                                     NTL::pow(P_induce,NTL::to_RR(n_to_induce)) *
                                     NTL::pow(NTL::RR(1)-P_induce,NTL::to_RR(n-t-n_to_induce));
             int n_to_correct = (int)t + n_to_induce - tau;
             NTL::RR prob_correct_n = NTL::to_RR(binomial_wrapper(t,n_to_correct));
                     prob_correct_n *= NTL::pow(P_correct,NTL::to_RR(n_to_correct));

                     prob_correct_n *= NTL::pow(NTL::RR(1)-P_correct,NTL::to_RR((int)t-n_to_correct)); /*unsigned exp?*/
             DFR += prob_correct_n*prob_induce_n;
         }
       }
       DFR = NTL::RR(1) - DFR;
       if (DFR < best_DFR){
          best_DFR = DFR;
          best_threshold = b;
       }
    }
    return std::make_pair(best_DFR,best_threshold);
}

// find minimum p which, asymptotically, corrects all errors
// search performed via binary search as the DFR is decreasing monot.
// in of p
uint64_t Findpth(const uint64_t n_0,
                 const uint64_t d_v_prime,
                 const uint64_t t){

    unsigned int prime_idx = 0, prime_idx_prec;
    uint64_t p = proper_primes[prime_idx];
    while(p < d_v_prime || p < t ){
          prime_idx++;
          p=proper_primes[prime_idx];
    }

    uint64_t hi, lo;
    lo = prime_idx;
    hi = PRIMES_NO;
    prime_idx_prec = lo;

    uint64_t limit_error_num = t;
    while(hi-lo > 1){
        prime_idx_prec = prime_idx;
        prime_idx = (lo+hi)/2;
        p = proper_primes[prime_idx];
        // compute number of remaining errors after +infty iters
        limit_error_num = t;
        uint64_t current_error_num;
//        std::cout << "using p:"<< p << ", errors dropping as ";
        do {
            current_error_num = limit_error_num;
            limit_error_num = FindNextNumErrors(n_0, p, d_v_prime, current_error_num);
//           std::cout << limit_error_num << " ";
           } while ( 
                     (limit_error_num != current_error_num) && 
                     (limit_error_num != 0)
                   );
//        std::cout << std::endl;
        if (limit_error_num > 0){
            lo = prime_idx;
        } else {
            hi = prime_idx;
        }
    }
    if(limit_error_num == 0) {
        return proper_primes[prime_idx];
    }
    return proper_primes[prime_idx_prec];
}
