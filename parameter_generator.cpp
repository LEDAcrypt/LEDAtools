#include <NTL/ZZ.h>
#include <cmath>
#include <cstdint>
#include <vector>

#define NUM_BITS_REAL_MANTISSA 128
#define IGNORE_DECODING_COST 0
#define SKIP_BJMM 1
#define SKIP_MMT 1
#define LOG_COST_CRITERION 1

#include "binomials.hpp"
#include "bit_error_probabilities.hpp"
#include "isd_cost_estimate.hpp"
#include "logging.hpp"
#include "partitions_permanents.hpp"
#include "proper_primes.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>

uint32_t estimate_t_val(const uint32_t c_sec_level, const uint32_t q_sec_level,
                        const uint32_t n_0, const uint32_t p) {
  double achieved_c_sec_level = c_sec_level;
  double achieved_q_sec_level = q_sec_level;
  uint32_t lo = 1, t, t_prec;
  uint32_t hi;
  hi = p < 4 * c_sec_level ? p : 4 * c_sec_level;
  t = lo;
  t_prec = lo;
  while (hi - lo > 1) {
    t_prec = t;
    t = (lo + hi) / 2;
    std::cerr << "testing t " << t << std::endl;
    achieved_c_sec_level =
        c_isd_log_cost(n_0 * p, ((n_0 - 1) * p), t, p, 0, true).value;
    achieved_q_sec_level =
        q_isd_log_cost(n_0 * p, ((n_0 - 1) * p), t, p, 0, true).value;
    if ((achieved_c_sec_level >= c_sec_level) &&
        (achieved_q_sec_level >= q_sec_level)) {
      hi = t;
    } else {
      lo = t;
    }
  }
  if ((achieved_c_sec_level >= c_sec_level) &&
      (achieved_q_sec_level >= q_sec_level)) {
    return t;
  }
  return t_prec;
}

int ComputeDvMPartition(const uint64_t d_v_prime, const uint64_t n_0,
                        std::vector<uint64_t> &mpartition,
                        uint64_t &d_v) {
  d_v = floor(sqrt(d_v_prime));
  d_v = (d_v & 0x01) ? d_v : d_v + 1;
  uint64_t m = ceil((double)d_v_prime / (double)d_v);

  int partition_ok;
  partition_ok = FindmPartition(m, mpartition, n_0);

  while (!partition_ok && (d_v_prime / d_v) >= n_0) {
    d_v += 2;
    m = ceil((double)d_v_prime / (double)d_v);
    partition_ok = FindmPartition(m, mpartition, n_0);
  }
  return partition_ok;
}

uint64_t
estimate_dv(const uint32_t c_sec_level, // expressed as
            const uint32_t q_sec_level, const uint32_t n_0, const uint32_t p,
            std::vector<uint64_t> &mpartition) {
  double achieved_c_sec_level = 0.0;
  double achieved_q_sec_level = 0.0;
  double achieved_c_enum_sec_level = 0.0;
  double achieved_q_enum_sec_level = 0.0;

  NTL::ZZ keyspace;

  uint32_t lo = 1, d_v_prime, hi;
  uint64_t d_v, d_v_prec = 0;
  int found_dv_mpartition = 0;
  // recalling that the weight of the sought codeword in a KRA is
  // d_c_prime = n_0 * d_v_prime, d_c_prime < p
  // d_v_prime should not be greater than p/n_0
  hi = (p / n_0) < 4 * c_sec_level ? (p / n_0) : 4 * c_sec_level;
  d_v_prime = lo;
  d_v = (int)sqrt(lo);

  while (hi - lo > 1) {
    d_v_prec = d_v;
    d_v_prime = (lo + hi) / 2;
    found_dv_mpartition = ComputeDvMPartition(d_v_prime, n_0, mpartition, d_v);
    if (found_dv_mpartition) {
      keyspace = 1;
      for (int i = 0; i < (int)n_0; i++) {
        keyspace *= binomial_wrapper(p, mpartition[i]);
      }
      keyspace = NTL::power(keyspace, n_0);
      keyspace += NTL::power(binomial_wrapper(p, d_v), n_0);
      achieved_c_enum_sec_level =
          NTL::conv<double>(log2_RR(NTL::to_RR(keyspace)));
      achieved_q_enum_sec_level = achieved_c_enum_sec_level / 2;
      if ((achieved_q_enum_sec_level >= q_sec_level) &&
          (achieved_c_enum_sec_level >= c_sec_level)) {
        /* last parameter indicates a KRA, reduce margin by p due to
        quasi cyclicity */
        achieved_c_sec_level =
            c_isd_log_cost(n_0 * p, p, n_0 * d_v_prime, p, 1, true).value;
        achieved_q_sec_level =
            q_isd_log_cost(n_0 * p, p, n_0 * d_v_prime, p, 1, true).value;
      }
    }

    if ((found_dv_mpartition) && (achieved_q_enum_sec_level >= q_sec_level) &&
        (achieved_c_enum_sec_level >= c_sec_level) &&
        (achieved_c_sec_level >= c_sec_level) &&
        (achieved_q_sec_level >= q_sec_level)) {
      hi = d_v_prime;
    } else {
      lo = d_v_prime;
    }
  }

  if ((found_dv_mpartition) && (achieved_q_enum_sec_level >= q_sec_level) &&
      (achieved_c_enum_sec_level >= c_sec_level) &&
      (achieved_c_sec_level >= c_sec_level) &&
      (achieved_q_sec_level >= q_sec_level)) {
    return d_v;
  }
  return d_v_prec;
}

int main(int argc, char *argv[]) {

  if (argc != 6) {
    std::cout << "Code Parameter Computer for LEDA[kem|pkc]" << std::endl
              << " Usage " << argv[0]
              << " security_level_classic security_level_pq n_0 epsilon "
                 "starting_prime_lb"
              << std::endl;
    return -1;
  }
  uint32_t c_sec_level = atoi(argv[1]);
  uint32_t q_sec_level = atoi(argv[2]);
  uint32_t n_0 = atoi(argv[3]);
  float epsilon = atof(argv[4]);
  uint32_t starting_prime_lower_bound = atoi(argv[5]);

  std::cerr << "Computing the parameter set for security level classic:2^"
            << c_sec_level << " post-q:2^" << q_sec_level << " n_0 " << n_0
            << " epsilon " << epsilon << std::endl;

  uint64_t p, p_th, t, d_v_prime, d_v;
  std::vector<uint64_t> mpartition(n_0, 0);

  int current_prime_pos = 0;
  while (proper_primes[current_prime_pos] < starting_prime_lower_bound) {
    current_prime_pos++;
  }
  p_th = proper_primes[current_prime_pos];

  InitBinomials();
  NTL::RR::SetPrecision(NUM_BITS_REAL_MANTISSA);
  pi = NTL::ComputePi_RR();

  /* since some values of p may yield no acceptable partitions for m, binary
   * search on p is not feasible. Fall back to fast increase of the value of
   * p exploiting the value of the p expected to be correcting the required t
   * errors */

  std::cout << "finding parameters" << std::endl;
  do {
    /* estimate the current prime as the closest
     * to the previous p_th * (1+epsilon) */
    uint32_t next_prime = ceil(p_th * (1.0 + epsilon));
    current_prime_pos = 0;
    while (proper_primes[current_prime_pos] < next_prime) {
      current_prime_pos++;
    }
    p = proper_primes[current_prime_pos];
    std::cout << " -- testing p: " << p << std::endl;

    // Estimate number of errors to ward off ISD decoding
    t = estimate_t_val(c_sec_level, q_sec_level, n_0, p);
    std::cout << " -- found t: " << t << std::endl;

    /* Estimate H*Q density to avoid key recovery via ISD and enumeration
     * of H and Q */
    d_v = estimate_dv(c_sec_level, q_sec_level, n_0, p, mpartition);
    std::cout << " -- found d_v: " << d_v << std::endl;

    // Estimate the bit flipping thresholds and correction capability
    d_v_prime = 0;
    for (int i = 0; i < (int)n_0; i++) {
      d_v_prime += mpartition[i];
    }
    d_v_prime = d_v * d_v_prime;
    p_th = Findpth(n_0, d_v_prime, t);
    std::cout << " -- p should be at least " << (1.0 + epsilon) * p_th
              << "to correct the errors" << std::endl;
  } while ((p <= (1.0 + epsilon) * p_th) && (current_prime_pos < PRIMES_NO));

  std::cout << "refining parameters" << std::endl;

  std::optional<uint32_t> p_ok;
  uint64_t t_ok, d_v_ok;
  std::vector<uint64_t> mpartition_ok(n_0, 0);
  /* refinement step taking into account possible invalid m partitions */

  do {
    p = proper_primes[current_prime_pos];
    std::cout << " -- testing p: " << p << std::endl;

    // Estimate number of errors to ward off ISD decoding
    t = estimate_t_val(c_sec_level, q_sec_level, n_0, p);
    std::cout << " -- found t: " << t << std::endl;

    /* Estimate H*Q density to avoid key recovery via ISD and enumeration
     * of H and Q */
    d_v = estimate_dv(c_sec_level, q_sec_level, n_0, p, mpartition);
    std::cout << " -- found d_v: " << d_v << std::endl;

    // Estimate the bit flipping thresholds and correction capability
    d_v_prime = 0;
    for (int i = 0; i < (int)n_0; i++) {
      d_v_prime += mpartition[i];
    }
    d_v_prime = d_v * d_v_prime;
    p_th = Findpth(n_0, d_v_prime, t);
    std::cout << " -- the threshold value for p to be correcting errors is "
              << p_th << std::endl;

    if (p > (1.0 + epsilon) * p_th) { // store last valid parameter set
      std::cout << " -- p is at least " << (1.0 + epsilon) * p_th
                << "; it corrects the errors" << std::endl;
      p_ok = p;
      t_ok = t;
      d_v_ok = d_v;
      for (unsigned i = 0; i < n_0; i++) {
        mpartition_ok[i] = mpartition[i];
      }
    }
    current_prime_pos--;
  } while ((p > (1.0 + epsilon) * p_th) && (current_prime_pos > 0));

  if (!p_ok || !d_v_ok) {
    spdlog::error("Error: One or more variables are not initialized.");
    throw std::runtime_error("One or more variables are not initialized.");
  } else {
    spdlog::info("parameter set found: p={}, t={}, d_v={}, mpartition={}",
                 optional_to_string(p_ok), t_ok, optional_to_string(d_v_ok),
                 array_to_string(mpartition_ok));
  }
    //         std::cout
    //     << " p:" << p_ok << " t: " << t_ok;
    // std::cout << " d_v : " << d_v_ok << " mpartition: [ ";
    // for (unsigned i = 0; i < n_0; i++) {
    //   std::cout << mpartition_ok[i] << " ";
    // }
    // std::cout << " ]" << std::endl;
    return 0;
  }
