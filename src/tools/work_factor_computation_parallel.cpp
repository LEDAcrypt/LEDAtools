#include <NTL/ZZ.h>
#include <cstdint>
#include <iomanip> // For std::setprecision
#include <iostream>
#include <nlohmann/json.hpp>
#include <omp.h>
#include <string>
#include <binomials.hpp>
#include <ctime>
#include <fstream>
#include <iostream>
#include <isd_cost_estimate.hpp>
// #include <logging.hpp>
#include <unordered_set>
#include <fmt/core.h>

#include "globals.hpp"

#define NUM_BITS_REAL_MANTISSA 1024
#define IGNORE_DECODING_COST 0
// #define EXPLORE_REPRS

void to_json(nlohmann::json &j, const Result &r) {
  j = nlohmann::json{
    {"alg_name", r.alg_name}, {"params", r.params}, {"value", r.value}, {"gje_cost", r.gje_cost}, {"list_size", r.list_size}};
}

void from_json(const nlohmann::json &j, Result &r) {
  j.at("alg_name").get_to(r.alg_name);
  j.at("params").get_to(r.params);
  j.at("value").get_to(r.value);
  j.at("gje_cost").get_to(r.gje_cost);
}

int main() {

  // Logger::LoggerManager::getInstance().setup_logger(
  //     "binomials", spdlog::level::info, spdlog::level::debug);
  // Logger::LoggerManager::getInstance().setup_logger(
  //     "isd_cost_estimate", spdlog::level::info, spdlog::level::debug);

  const std::string input_isd_values = "out/isd_values.json";
  std::ifstream file(input_isd_values);

  // Check if the file is open 
  if (!file.is_open()) {
    std::cerr << "Could not open the input file " << input_isd_values  << std::endl;
    return 1;
  }

  // Parse the JSON content
  nlohmann::json j;
  file >> j;

  NTL::RR::SetPrecision(NUM_BITS_REAL_MANTISSA);
  
  InitBinomials();
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
        OUT_DIR_RESULTS + fmt::format("{:06}_{:06}_{:03}.json", n, r, t);

    nlohmann::json out_values;

    Result current_c_res;
    Result current_q_res;

    for (bool is_kra : is_kra_values) {
      // spdlog::info("Processing n {}, k {}, t {}, qc_block_size {}, is_kra {}, "
      //              "is_red_factor_applied {}",
      //              n, k, t, qc_block_size, is_kra, is_red_factor_applied);
      current_c_res =
          c_isd_log_cost(n, k, t, qc_block_size, is_kra, is_red_factor_applied,
                         std::unordered_set<Algorithm>{Stern});
      current_q_res =
          q_isd_log_cost(n, k, t, qc_block_size, is_kra, is_red_factor_applied,
                         std::unordered_set<QuantumAlgorithm>{Q_Lee_Brickell});
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
