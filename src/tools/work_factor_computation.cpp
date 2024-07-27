#include <NTL/ZZ.h>
#include <cstdint>

#define NUM_BITS_REAL_MANTISSA 1024
#define IGNORE_DECODING_COST 0
// #define EXPLORE_REPRS

#include "binomials.hpp"
#include "isd_cost_estimate.hpp"
#include <cstdlib>
#include <iostream>

int main(int argc, char *argv[]) {
  if (argc != 7) {
    std::cout
        << "Work factor computation for ISD" << std::endl
        << " Usage " << argv[0]
        << " <codeword_size> <code_dimension> <number_of_errors> "
           "<qc_block_size> <is_kra>"
        << std::endl
        << "<qc_block_size> = 1 implies a non QC code " << std::endl
        << "<is_kra> = the attack is a key recovery attack on a QC-[L|M]DPC "
        << std::endl
        << "<is_red_factor_applied> = if the quasi-cyclic reduction factor "
           "cost should be applied"
        << std::endl;
    return -1;
  }

  /* reduce by a factor matching the QC block size */

  
  InitBinomials();
  NTL::RR::SetPrecision(NUM_BITS_REAL_MANTISSA);
  pi = NTL::ComputePi_RR();
  uint32_t n = atoi(argv[1]);
  uint32_t k = atoi(argv[2]);
  uint32_t t = atoi(argv[3]);
  uint32_t qc_block_size = atoi(argv[4]);
  uint32_t is_kra = atoi(argv[5]);
  uint32_t is_red_factor_applied = atoi(argv[6]);

  std::cout << " Input params: " << std::endl
            << "- <codeword_size>: " << n << std::endl
            << "- <code_dimension>: " << k << std::endl
            << "- <number of errors>: " << t << std::endl
            << "- <qc block size>: " << qc_block_size << std::endl
            << "- <is_kra>: " << is_kra << std::endl
            << "- <is_red_factor_applied>: " << is_red_factor_applied << std::endl;

  std::cout << "Minimum classic cost :"
            << c_isd_log_cost(
                   n, k, t, qc_block_size, is_kra, is_red_factor_applied,
                   std::unordered_set<Algorithm>{Prange, Lee_Brickell, Leon,
                                                 Stern, Finiasz_Sendrier, MMT,
                                                 BJMM})
                   .value
            << " Minimum quantum cost :"
            << q_isd_log_cost(n, k, t, qc_block_size, is_kra,
                              is_red_factor_applied,
                              std::unordered_set<QuantumAlgorithm>{
                                  Q_Lee_Brickell, Q_Stern})
                   .value;
  if (is_red_factor_applied && qc_block_size != 1)
    std::cout << " (including qc_effects) ";
  std::cout << std::endl;
  return 0;
}
