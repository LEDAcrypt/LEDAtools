#include <NTL/ZZ.h>
#include <cstdint>
#include <fstream>
#include <iomanip> // For std::setprecision
#include <iostream>
#include <nlohmann/json.hpp>
#include <omp.h>
#include <sstream>
#include <string>

#include "binomials.hpp"
#include "isd_cost_estimate.hpp"
#include "logging.hpp"
#include "globals.hpp"
#include <iostream>
#include <format> // Requires C++20

#define NUM_BITS_REAL_MANTISSA 1024
#define IGNORE_DECODING_COST 0
// #define EXPLORE_REPRS

void to_json(nlohmann::json &j, const Result &r) {
  j = nlohmann::json{
      {"alg_name", r.alg_name}, {"params", r.params}, {"value", r.value}};
}

void from_json(const nlohmann::json &j, Result &r) {
  j.at("alg_name").get_to(r.alg_name);
  j.at("params").get_to(r.params);
  j.at("value").get_to(r.value);
}

int main() {
  // Configure the logger

  configure_logger(std::nullopt);

  // TODO take from input?
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
  std::filesystem::path dirPath(OUT_DIR_RESULTS);
  // Check if the directory exists
  if (!std::filesystem::exists(dirPath)) {
    // Try to create the directory, including parent directories
    if (std::filesystem::create_directories(dirPath)) {
      std::cout << "Directory created successfully: " << OUT_DIR_RESULTS
                << std::endl;
    } else {
      std::cerr << "Failed to create directory: " << OUT_DIR_RESULTS
                << std::endl;
      return 1; // Return an error code
    }
  }
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

    std::string filename =
        OUT_DIR_RESULTS + fmt::format("/{:06}_{:06}_{:03}.json", n, r, t);

    nlohmann::json out_values;

    Result current_c_res;
    Result current_q_res;

    for (bool is_kra : is_kra_values) {
      spdlog::info("Processing n {}, k {}, t {}, qc_block_size {}, is_kra {}, "
                   "is_red_factor_applied {}",
                   n, k, t, qc_block_size, is_kra, is_red_factor_applied);
      current_c_res =
          c_isd_log_cost(n, k, t, qc_block_size, is_kra, is_red_factor_applied);
      current_q_res =
          q_isd_log_cost(n, k, t, qc_block_size, is_kra, is_red_factor_applied);
      std::string is_kra_name = is_kra ? "KRA": "MRA";
      out_values[is_kra_name]["C"] = current_c_res;
      out_values[is_kra_name]["Q"] = current_q_res;
    }

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
  return 0;
}
