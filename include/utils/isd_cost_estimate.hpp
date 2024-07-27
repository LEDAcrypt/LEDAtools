#pragma once
#include <NTL/RR.h>
#include <NTL/ZZ.h>
#include <map>
#include <optional>
#include <string>
#include <unordered_set>

struct Result {
  std::string alg_name;
  std::map<std::string, int> params;
  double value;
  double gje_cost;
  double list_size;
};

// Plain does not apply qc reductions
enum class QCAttackType { KRA1, KRA2, KRA3, MRA, Plain, Count};

enum class Algorithm {
  Prange,
  Lee_Brickell,
  Leon,
  Stern,
  Finiasz_Sendrier,
  MMT,
  BJMM,
  // Add more algorithms here
  Count,
};
enum class QuantumAlgorithm {
  Q_Lee_Brickell,
  Q_Stern, // NOTE no circuit available
};

/***************************Classic ISDs***************************************/

const NTL::RR log_probability_k_by_k_is_inv(const NTL::RR &k);
const NTL::RR probability_k_by_k_is_inv(const NTL::RR &k);
const NTL::RR classic_rref_red_cost(const NTL::RR &n, const NTL::RR &r);

// Classic
Result c_isd_log_cost(const uint32_t n, const uint32_t k, const uint32_t t,
                      const uint32_t qc_order, QCAttackType attack,
                      const bool compute_qc_reduction_factor,
                      std::unordered_set<Algorithm> algs);

double get_qc_red_factor_log(const uint32_t qc_order, const uint32_t n0,
                             QCAttackType attack);

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
                      const uint32_t qc_order, QCAttackType attack,
                      const bool compute_qc_reduction_factor,
                      std::unordered_set<QuantumAlgorithm> algs);

Result isd_log_cost_quantum_LB(const uint32_t n, const uint32_t k,
                               const uint32_t t);
Result isd_log_cost_quantum_Stern(const uint32_t n, const uint32_t k,
                                  const uint32_t t);
