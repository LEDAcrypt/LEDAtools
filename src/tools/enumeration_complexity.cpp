#include <NTL/ZZ.h>
#include <cstdint>
#include <cmath>

#define NUM_BITS_REAL_MANTISSA 54

#include "binomials.hpp"
#include <cmath>

int main(int argc, char* argv[]){
  if(argc < 5){
     std::cout << "Complexity estimation of key enumeration" << std::endl << " Usage " 
               << argv[0] << " <p> <d_v> <n_0> <m_0> <m_1> <m_2> <m_3>" << std::endl << 
               "All the non existing m_i should be passed as zeroes " << std::endl;
    return -1;
  }

  
  InitBinomials();
  NTL::RR::SetPrecision(NUM_BITS_REAL_MANTISSA);
  uint32_t p = atoi(argv[1]);
  uint32_t d_v = atoi(argv[2]);
  uint32_t n_0 = atoi(argv[3]);
  uint32_t m[4];
  m[0] = atoi(argv[4]);
  m[1] = atoi(argv[5]);
  m[2] = atoi(argv[6]);
  m[3] = atoi(argv[7]);

  NTL::RR Henum, Qenum;

  Henum = lnBinom(NTL::RR(p),NTL::RR(d_v));
  Henum = Henum*NTL::RR(n_0) / NTL::log(NTL::RR(2));

  Qenum = NTL::RR(n_0)* ( lnBinom(NTL::RR(p),NTL::RR(m[0])) +
                          lnBinom(NTL::RR(p),NTL::RR(m[1])) +
                          lnBinom(NTL::RR(p),NTL::RR(m[2])) +
                          lnBinom(NTL::RR(p),NTL::RR(m[3])) )/ NTL::log(NTL::RR(2));

  std::cout << "H enum classic/quantum cost :" << Henum << "  " << (Henum/NTL::RR(2)) << std::endl;
  std::cout << "Q enum classic/quantum cost :" << Qenum << "  " << (Qenum/NTL::RR(2)) << std::endl;
  return 0;
}
