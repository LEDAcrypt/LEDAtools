#pragma once
#include <cstdint>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/matrix.h>

/* binomials are precomputed up to MAX_N-choose-MAX_T */
#define MAX_N 2000
#define MAX_T 300

#define LOW_K_MAX_N 10000
#define LOW_K_MAX_T 10

const NTL::RR nat_log_2 = NTL::log(NTL::RR(2));

static inline NTL::RR log2_RR(NTL::RR v){
    return NTL::log(v)/nat_log_2;
}

NTL::Mat<NTL::ZZ> binomial_table; /*contains all binomials up to MAX_N-choose-MAX_T */
NTL::Mat<NTL::ZZ> low_k_binomial_table; /*contains all binomials up to LOW_K_MAX_N-choose-LOW_K_MAX_T */

NTL::RR pi;
/*NOTE: NTL allows to access matrices as 1- based with Matlab notation */
void InitBinomials(){
    std::cerr << "Precomputing n-choose-t up to n: " << MAX_N << 
                                              " t: " << MAX_T << std::endl;
    binomial_table.SetDims(MAX_N+1,MAX_T+1);
    binomial_table[0][0] = NTL::ZZ(1);
    for (unsigned i = 1 ; i <= MAX_N; i++){
        binomial_table[i][0] = NTL::ZZ(1);
        binomial_table[i][1] = NTL::ZZ(i);
        for(unsigned j=2 ; (j <= i) && (j <= MAX_T) ; j++){
            binomial_table[i][j] = binomial_table[i][j-1] * NTL::ZZ(i-j+1) / NTL::ZZ(j);
        }
    }
    std::cerr << "Precomputing low n-choose-t up to n: " << LOW_K_MAX_N << 
                                                       " t: " << LOW_K_MAX_T << std::endl;
    low_k_binomial_table.SetDims(LOW_K_MAX_N+1,LOW_K_MAX_T+1);
    low_k_binomial_table[0][0] = NTL::ZZ(1);
    for (unsigned i = 0 ; i <= LOW_K_MAX_N; i++){
        low_k_binomial_table[i][0] = NTL::ZZ(1);
        low_k_binomial_table[i][1] = NTL::ZZ(i);
        for(unsigned j=2 ; (j <= i) && (j <= LOW_K_MAX_T) ; j++){
            low_k_binomial_table[i][j] = low_k_binomial_table[i][j-1] * NTL::ZZ(i-j+1) / NTL::ZZ(j);
        }
    }
    std::cerr << "done" << std::endl;
    pi = NTL::ComputePi_RR();
}



NTL::RR lnFactorial(NTL::RR n){
    /* log of Stirling series approximated to the fourth term 
     * n log(n) - n + 1/2 log(2 \pi n) + log(- 139/(51840 n^3) +
     * + 1/(288 n^2) + 1/(12 n) + 1) */
    return n * NTL::log(n) - n + 0.5 * NTL::log(2*pi*n) + 
           NTL::log( - NTL::RR(139)/(n*n*n * 51840) + 
           NTL::RR(1)/(n*n*288) + 
           NTL::RR(1)/(n*12) + 
           1);
}

NTL::RR lnBinom(NTL::RR n, NTL::RR k){
    if ( (k == NTL::RR(0) ) || (k == n) ) {
        return NTL::RR(0);
    }
    return lnFactorial(n) - (lnFactorial(k) + lnFactorial(n-k) );
}


NTL::ZZ binomial_wrapper(long n, long k){
    if(k>n) return NTL::ZZ(0);
    /* employ memoized if available */
    if ((n <= MAX_N) && (k < MAX_T)){
        return binomial_table[n][k];
    }
    if ((n <= LOW_K_MAX_N) && (k < LOW_K_MAX_T)){
        return low_k_binomial_table[n][k];
    }
    
    /* shortcut computation for fast cases (k < 10) where 
     * Stirling may not provide good approximations */
    if (k < 10) {
        NTL::ZZ result = NTL::ZZ(1);
        for(int i = 1 ; i <= k; i++){
            result = (result * (n+1-i))/i;
        }
        return result;
    }
    /*Fall back to Stirling*/
    return NTL::conv<NTL::ZZ>( NTL::exp( lnBinom(NTL::RR(n),NTL::RR(k)) ));
}
