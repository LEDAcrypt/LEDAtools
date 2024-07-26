#define NUM_BITS_REAL_MANTISSA 128
#include <cstdint>
#include <cmath>
#include <NTL/ZZ.h>

#include "utils/binomials.hpp"

int main(int argc, char* argv[]){
 if(argc != 3){
     std::cout << "Calculator to derive the length of the encodable bit string via CW-enc" << std::endl << " Usage " 
               << argv[0] << " <codeword_size> <number_of_errors> " << std::endl;
    return -1;
 }

 InitConstants();
 InitBinomials();
 NTL::RR::SetPrecision(NUM_BITS_REAL_MANTISSA);
 uint32_t n = atoi(argv[1]);
 uint32_t t = atoi(argv[2]);
 /* reduce by a factor matching the QC block size */
 NTL::RR encodable_length;
 encodable_length = lnBinom(NTL::to_RR(n), NTL::to_RR(t))/NTL::log(NTL::RR(2));

 NTL::RR d = NTL::to_RR( 0.69315 * ((double)n - ( (double)t - 1.0)/2.0) /((double) t) );

 std::cout << "Maximum safely encoded: " << t*NTL::conv<unsigned long int>(NTL::floor(NTL::log(d)/NTL::log(NTL::RR(2))+1)) << std::endl;
 std::cout << "#define MAX_ENCODABLE_BIT_SIZE_CW_ENCODING (" << NTL::conv<unsigned long int>(encodable_length) << ")" ;
  std::cout << std::endl;
  return 0;
}
