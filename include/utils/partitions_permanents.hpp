#include <cstdint>
#include <iostream>
#include <vector>

/**
 * @brief Computes the permanent of a circulant matrix.
 *
 * @param mpartition Array of integers representing the partition.
 * @param n_0 Size of the partition.
 * @return uint64_t Permanent of the circulant matrix.
 */
uint64_t ComputePermanent(int64_t mpartition[], const uint64_t n_0);

/**
 * @brief Finds a partition of m with length n_0.
 *
 * @param m The integer to partition.
 * @param mpartition Vector to store the resulting partition.
 * @param n_0 Length of the partition.
 * @return int 1 if a good partition is found, 0 otherwise.
 */
int FindmPartition(const uint64_t m, std::vector<uint64_t> &mpartition,
                   const uint64_t n_0);
