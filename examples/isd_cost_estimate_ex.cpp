#include <cstdint> // for uint32_t
#include <cstdlib>
#include <iostream>
#include <isd_cost_estimate.hpp>
#include <map>
#include <ostream>
#include <spdlog/spdlog.h>
#include <string>
#include <vector>

// #include "logging.hpp"

struct Cost {
  std::string algorithm;
  std::string type; // CFP1, CFP2, CFP3, SDP
  bool is_quantum;
  double time_complexity;
  double space_complexity;
};

struct Value {
  uint32_t codeword_size;
  uint32_t code_dimension;
  uint32_t number_of_errors;
  uint32_t qc_block_size;
  bool is_kra;
  bool is_red_fac;
  std::map<std::string, Cost> costs;
};

void displayCost(const Cost &cost) {
  std::cout << "  Algorithm: " << cost.algorithm << '\n';
  std::cout << "  Type: " << cost.type << '\n';
  std::cout << "  Is Quantum: " << (cost.is_quantum ? "Yes" : "No") << '\n';
  std::cout << "  Time Complexity: " << cost.time_complexity << '\n';
  std::cout << "  Space Complexity: " << cost.space_complexity << '\n';
}

// Function to display a Value object
void displayValue(const Value &value) {
  std::cout << "Value:\n";
  std::cout << "  Codeword Size: " << value.codeword_size << '\n';
  std::cout << "  Code Dimension: " << value.code_dimension << '\n';
  std::cout << "  Number of Errors: " << value.number_of_errors << '\n';
  std::cout << "  QC Block Size: " << value.qc_block_size << '\n';
  std::cout << "  Is KRA: " << (value.is_kra ? "Yes" : "No") << '\n';
  std::cout << "  Is Reduction factor applied: "
            << (value.is_red_fac ? "Yes" : "No") << '\n';

  std::cout << "Costs:\n";
  for (const auto &costPair : value.costs) {
    displayCost(costPair.second);
    std::cout << "-----\n";
  }
}

int main() {
  std::cout<< "logger setted up" << std::endl;

  std::vector<Value> values;
  Cost pra = {"Prange", "", false, 171.3, 0.0};
  Cost lbr = {"Lee-Brickell", "", false, 158.4, 0.0};
  Cost leo = {"Leon", "", false, 154.4, 0.0};
  Cost ste = {"Stern", "", false, 147.4, 0.0};
  Cost fis = {"Fin-Send", "", false, 147.4, 0.0};

  std::map<std::string, Cost> costs = {{"Prange", pra},
                                       {"Lee-Brickell", lbr},
                                       {"Leon", leo},
                                       {"Stern", ste},
                                       {"Fin-Send", fis}};

  Value val = {24646, 12323, 142, 12323, true, true, costs};

  // displayValue(val);

  Result current_res, min_res;
  uint32_t n = val.codeword_size;
  uint32_t k = val.code_dimension;
  uint32_t t = val.number_of_errors;
  double qc_red_fac = get_qc_red_factor_log(val.qc_block_size, n-k, QCAttackType::KRA3);
  double diff;
  std::string name;
  std::cout << "qc_red_fac " << qc_red_fac << std::endl;
  for (const auto &algo : std::unordered_set<Algorithm>{
           Algorithm::Prange, Algorithm::Lee_Brickell, Algorithm::Leon,
           Algorithm::Stern, Algorithm::Finiasz_Sendrier}) {
    switch (algo) {
    case Algorithm::Prange:
      current_res = isd_log_cost_classic_Prange(n, k, t);
      name = "Prange";
      break;
    case Algorithm::Lee_Brickell:
      current_res = isd_log_cost_classic_LB(n, k, t);
      name = "Lee-Brickell";
      break;
    case Algorithm::Leon:
      current_res = isd_log_cost_classic_Leon(n, k, t);
      name = "Leon";
      break;
    case Algorithm::Stern:
      current_res = isd_log_cost_classic_Stern(n, k, t);
      name = "Stern";
      break;
    case Algorithm::Finiasz_Sendrier:
      current_res = isd_log_cost_classic_FS(n, k, t);
      name = "Fin-Send";
      break;
    case Algorithm::MMT:
      current_res = isd_log_cost_classic_MMT(n, k, t);
      name = "MMT ";
      break;
    case Algorithm::BJMM:
      current_res = isd_log_cost_classic_BJMM(n, k, t);
      name = "BJMM ";
      break;
    default:
      std::cerr << "Unknown algorithm\n";
      break;
    }
    current_res.value -= qc_red_fac;
    diff = std::abs(costs[name].time_complexity - current_res.value);
    std::cout << name << ". Obtained: " << current_res.value
              << " Expected: " << costs[name].time_complexity
              << " Diff: " << diff << std::endl;
    if (diff >= 1.0) {
      std::cerr << "WARNING: huge diff";
    }
  }

  return 0;
}
