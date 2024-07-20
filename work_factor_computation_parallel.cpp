#include <NTL/ZZ.h>
#include <cstdint>
#include <fstream>
#include <iomanip> // For std::setprecision
#include <iostream>
#include <nlohmann/json.hpp>
#include <omp.h>
#include <sstream>
#include <string>

#define NUM_BITS_REAL_MANTISSA 1024
#define IGNORE_DECODING_COST 0
// #define EXPLORE_REPRS

#include "binomials.hpp"
#include "isd_cost_estimate.hpp"
#include <iostream>


int main(int argc, char *argv[]) {
  std::ifstream file("out/isd_values.json");

  // Check if the file is open
  if (!file.is_open()) {
    std::cerr << "Could not open the file!" << std::endl;
    return 1;
  }

  // Parse the JSON content
  nlohmann::json j;
  file >> j;

  InitBinomials();
  NTL::RR::SetPrecision(NUM_BITS_REAL_MANTISSA);
  pi = NTL::ComputePi_RR();
  bool is_kra_values[] = {true, false};
  // Iterate over the list of entries
  #pragma omp parallel for
  for (const auto &entry : j) {
    uint32_t n = entry["n"];
    uint32_t r = entry["r"];
    uint32_t k = n - r;
    uint32_t t = entry["t"];
    uint32_t qc_block_size = entry["prime"];
    bool is_red_factor_applied = true;
    // int n0 = entry["n0"];
    // int v = entry["v"];
    // int lambd = entry["lambd"];

    // Output the data
    // std::cout << "n: " << n << ", r: " << r << ", t: " << t << std::endl;

    for (bool is_kra : is_kra_values) {
      double min_c_cost =
          c_isd_log_cost(n, k, t, qc_block_size, is_kra, is_red_factor_applied);
      double min_q_cost =
          q_isd_log_cost(n, k, t, qc_block_size, is_kra, is_red_factor_applied);
      nlohmann::json out_values;
      out_values["C2"] = min_c_cost;
      out_values["Q2"] = min_q_cost;

      std::ostringstream oss;
      oss << std::setw(6) << std::setfill('0') << n << "_" << std::setw(6)
          << std::setfill('0') << r << "_" << std::setw(3) << std::setfill('0')
          << t << "_" << std::setw(1) << is_kra;
      std::string filename = "out/" + oss.str() + ".json";

      // Write the JSON object to the file
      std::ofstream file(filename);
      if (file.is_open()) {
        file << std::fixed << std::setprecision(10)
             << out_values.dump(4); // Format JSON with indentation
        file.close();
        std::cout << "Data written to " << filename << std::endl;
      } else {
        std::cerr << "Could not open the file!" << std::endl;
      }
    }
  }
  return 0;
}
