#define NUM_BITS_REAL_MANTISSA 128
#include <cstdint>
#include <cmath>
#include <NTL/ZZ.h>

#include "binomials.hpp"

int main(int argc, char* argv[]){
 if(argc != 3){
     std::cout << "Calculator to derive the length of the encodable bit string via CW-enc" << std::endl << " Usage " 
               << argv[0] << " <codeword_size> <number_of_errors> " << std::endl;
    return -1;
 }

 InitBinomials();
 NTL::RR::SetPrecision(NUM_BITS_REAL_MANTISSA);
 pi = NTL::ComputePi_RR();
 uint32_t n = atoi(argv[1]);
 uint32_t t = atoi(argv[2]);
 /* reduce by a factor matching the QC block size */
 NTL::RR encodable_length;
 encodable_length = lnBinom(NTL::to_RR(n), NTL::to_RR(t));
 std::cout << "#define MAX_ENCODABLE_BIT_SIZE_CW_ENCODING (" << NTL::conv<unsigned long int>(encodable_length/NTL::log(NTL::RR(2))) << ")" ;
  std::cout << std::endl;
  return 0;
}
