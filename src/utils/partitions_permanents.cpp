#include <cstdint>

#include "partitions_permanents.hpp"

/*
 Permanent formulas for circulant matrices as obtained via Sage (macsyma)

+           +
| [m_0 m_1] |
| [m_1 m_0] | = (m_0 + m_1)^2 - 2*m_0*m_1

         	
+               +
| [m_0 m_1 m_2] |
| [m_2 m_0 m_1] |
| [m_1 m_2 m_0] | = (m_0 + m_1 + m_2)^3 - 3*(m_0 + m_1)*(m_0 + m_2)*(m_1 + m_2) +
                  3*m_0*m_1*m_2

+                   +
| [m_0 m_1 m_2 m_3] |
| [m_3 m_0 m_1 m_2] |
| [m_2 m_3 m_0 m_1] |
| [m_1 m_2 m_3 m_0] | =  (m_0 + m_1 + m_2 + m_3)^4 - 
                      4*(m_0 + m_1 + m_2)*(m_0 + m_1 + m_3)*(m_0 +m_2 + m_3)* 
                      (m_1 + m_2 + m_3) + 
                      2*(m_0 + m_2)^2*(m_1 + m_3)^2 + 4*(m_0 + m_1)*(m_0 + m_3)*
                      (m_1 + m_2)*(m_2 + m_3) - 4*m_0*m_1*m_2*m_3
*/

uint64_t ComputePermanent(int64_t mpartition[], 
                          const uint64_t n_0){
            uint64_t permanent = 0;                  
switch(n_0){
    case 2: 
        permanent =  (mpartition[0] + mpartition[1]) *
                     (mpartition[0] + mpartition[1]) -
                     2*mpartition[0]*mpartition[1];
        break;
    case 3: 
        permanent =   (mpartition[0] + mpartition[1] + mpartition[2]) *
                      (mpartition[0] + mpartition[1] + mpartition[2]) *
                      (mpartition[0] + mpartition[1] + mpartition[2]) -
                    3*(mpartition[0] + mpartition[1]) *
                      (mpartition[0] + mpartition[2]) *
                      (mpartition[1] + mpartition[2]) +
                    3* mpartition[0] *
                       mpartition[1] *
                       mpartition[2];
        break;
    case 4:
        permanent =  (mpartition[0] + mpartition[1] + mpartition[2] + mpartition[3]) *
                     (mpartition[0] + mpartition[1] + mpartition[2] + mpartition[3]) *
                     (mpartition[0] + mpartition[1] + mpartition[2] + mpartition[3]) *
                     (mpartition[0] + mpartition[1] + mpartition[2] + mpartition[3]) -
                   4*(mpartition[0] + mpartition[1] + mpartition[2]) *
                     (mpartition[0] + mpartition[1] + mpartition[3]) *
                     (mpartition[0] + mpartition[2] + mpartition[3]) *
                     (mpartition[1] + mpartition[2] + mpartition[3]) +
                   2*((mpartition[0] + mpartition[2])*(mpartition[0] + mpartition[2])) * 
                     ((mpartition[1] + mpartition[3])*(mpartition[1] + mpartition[3])) + 
                   4*(mpartition[0] + mpartition[1]) * 
                     (mpartition[0] + mpartition[3]) *
                     (mpartition[1] + mpartition[2]) *
                     (mpartition[2] + mpartition[3]) - 
                   4* mpartition[0] * 
                      mpartition[1] *
                      mpartition[2] *
                      mpartition[3];
        break;
    default: std::cout << "permanent not supported" << std::endl;
  }
  return permanent;
}

int FindmPartition(const uint64_t m,
                   std::vector<uint64_t> &mpartition,
                   const uint64_t n_0) {
  // Enumerate partitions of m with length n_0,
  // according to TAOCP, Vol 4 Fascicle 3b, Algorithm H
  // PRE : m >= n_0 >= 2
  if ((m < n_0) || (n_0 < 2)) {
    return 0;
  }

  int64_t mpartition_selected[n_0];
  int found_good_partition = 0;
  int64_t mpartition_tmp[n_0 + 1];
  mpartition_tmp[0] = m - (n_0 - 1);
  for (unsigned i = 1; i < n_0; i++) {
    mpartition_tmp[i] = 1;
  }
  // theoretically, mpartition_tmp[n_0] = -1 according to knuth
  mpartition_tmp[n_0] = -1;
  do {
    // visit the partition
    if (ComputePermanent(mpartition_tmp, n_0) % 2 == 1U) {
      for (unsigned i = 0; i < n_0; i++) {
        mpartition_selected[i] = (uint64_t)mpartition_tmp[i];
      }
      found_good_partition = 1;
    }
    if (mpartition_tmp[1] < mpartition_tmp[0] - 1) {
      // step H3: easy but very common case
      mpartition_tmp[0]--;
      mpartition_tmp[1]++;
    } else {
      // step H4
      int j = 2;
      int64_t s = mpartition_tmp[0] + mpartition_tmp[1] - 1;
      while (mpartition_tmp[j] >= mpartition_tmp[0] - 1) {
        s = s + mpartition_tmp[j];
        j++;
      }
      if (j >= (int)n_0) { // completed enumeration of partitions
        for (unsigned i = 0; i < n_0; i++) {
          mpartition[i] = (uint64_t)mpartition_selected[i];
        }
        return found_good_partition;
      } else {
        uint64_t x = mpartition_tmp[j] + 1;
        mpartition_tmp[j] = x;
        j--;
        while (j > 0) {
          mpartition_tmp[j] = x;
          s = s - x;
          j--;
        }
        mpartition_tmp[0] = s;
      }
    }
  } while (1);
}
