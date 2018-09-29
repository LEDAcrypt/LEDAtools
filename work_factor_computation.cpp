#include <NTL/ZZ.h>
#include <cstdint>
#include <cmath>

#define NUM_BITS_REAL_MANTISSA 128
#define IGNORE_DECODING_COST 0

#include "binomials.hpp"
#include "isd_cost_estimate.hpp"
#include <cmath>

int main(int argc, char* argv[]){
  if(argc != 6){
     std::cout << "Work factor computation for ISD" << std::endl << " Usage " 
               << argv[0] << " <codeword_size> <code_dimension> <number_of_errors> <qc_block_size> <is_kra>" << std::endl << 
               "<qc_block_size> = 1 implies a non QC code " << std::endl << 
               "<is_kra> = the attack is a KRA on [L|M]DPC " << std::endl;
    return -1;
  }

  InitBinomials();
  NTL::RR::SetPrecision(NUM_BITS_REAL_MANTISSA);
  pi = NTL::ComputePi_RR();
  uint32_t n = atoi(argv[1]);
  uint32_t k = atoi(argv[2]);
  uint32_t t = atoi(argv[3]);
  uint32_t qc_block_size = atoi(argv[4]);
  uint32_t is_kra = atoi(argv[5]);
  
  /* reduce by a factor matching the QC block size */
  std::cout << "Minimum classic cost :" << c_isd_log_cost(n,k,t,qc_block_size,is_kra) << " Minimum quantum cost :" << q_isd_log_cost(n,k,t,qc_block_size,is_kra);
  if(qc_block_size !=1) std::cout << " (including qc_effects) ";
  std::cout << std::endl;
  return 0;
}
