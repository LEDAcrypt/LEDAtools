#pragma once
#include <NTL/RR.h>
#include <NTL/ZZ.h>
#include <map>
#include <optional>
#include <string>

#define SKIP_PRANGE 0
#define SKIP_LB 0
#define SKIP_LEON 0
#define SKIP_STERN 0
#define SKIP_FS 0
#define SKIP_BJMM 0
#define SKIP_MMT 0
#define SKIP_Q_LB 0
#define SKIP_Q_STERN 0

struct Result {
  std::string alg_name;
  std::map<std::string, int> params;
  double value;
  double gje_cost;
  double list_size;
};

/***************************Classic ISDs***************************************/

const NTL::RR log_probability_k_by_k_is_inv(const NTL::RR &k);
const NTL::RR probability_k_by_k_is_inv(const NTL::RR &k);
const NTL::RR classic_rref_red_cost(const NTL::RR &n, const NTL::RR &r);

// Classic

Result c_isd_log_cost(const uint32_t n, const uint32_t k, const uint32_t t,
                      const uint32_t qc_order, const uint32_t is_kra,
                      const bool compute_qc_reduction_factor);

Result isd_log_cost_classic_Prange(const uint32_t n, const uint32_t k,
                                   const uint32_t t);
Result isd_log_cost_classic_LB(const uint32_t n, const uint32_t k,
                               const uint32_t t);
Result isd_log_cost_classic_Leon(const uint32_t n, const uint32_t k,
                                 const uint32_t t);
Result isd_log_cost_classic_Stern(const uint32_t n, const uint32_t k,
                                  const uint32_t t);
Result isd_log_cost_classic_FS(const uint32_t n, const uint32_t k,
                               const uint32_t t);
Result isd_log_cost_classic_MMT(const uint32_t n, const uint32_t k,
                                const uint32_t t);
Result isd_log_cost_classic_BJMM_approx(const uint32_t n, const uint32_t k,
                                        const uint32_t t);
Result isd_log_cost_classic_BJMM(const uint32_t n, const uint32_t k,
                                 const uint32_t t);

// Quantum
Result q_isd_log_cost(const uint32_t n, const uint32_t k, const uint32_t t,
                      const uint32_t qc_order, const uint32_t is_kra,
                      const bool compute_qc_reduction_factor);

Result isd_log_cost_quantum_LB(const uint32_t n, const uint32_t k,
                               const uint32_t t);
Result isd_log_cost_quantum_Stern(const uint32_t n, const uint32_t k,
                                  const uint32_t t);
