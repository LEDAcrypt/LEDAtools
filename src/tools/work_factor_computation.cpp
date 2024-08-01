#include <NTL/ZZ.h>
#include <atomic>
#include <binomials.hpp>
#include <cstdint>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip> // For std::setprecision
#include <iostream>
#include <isd_cost_estimate.hpp>
#include <nlohmann/json.hpp>
#include <omp.h>
#include <string>
// #include <logging.hpp>
#include <fmt/core.h>
#include <unordered_set>

#include "globals.hpp"

#define NUM_BITS_REAL_MANTISSA 1024
#define IGNORE_DECODING_COST 0
// #define EXPLORE_REPRS

void to_json(nlohmann::json &j, const Result &r) {
  j = nlohmann::json{{"alg_name", r.alg_name},
                     {"params", r.params},
                     {"value", r.value},
                     {"gje_cost", r.gje_cost},
                     {"list_size", r.list_size}};
}

void from_json(const nlohmann::json &j, Result &r) {
  j.at("alg_name").get_to(r.alg_name);
  j.at("params").get_to(r.params);
  j.at("value").get_to(r.value);
  j.at("gje_cost").get_to(r.gje_cost);
}

int handle_plain(const std::string args) {
  std::istringstream argStream(args);
  std::string token;
  std::vector<int> values;
  while (std::getline(argStream, token, ',')) {
    values.push_back(std::stoi(token));
  }

  if (values.size() != 4) {
    std::cerr << "Expected 4 comma-separated values, but got " << values.size()
              << std::endl;
    return 1;
  }

  int n = values[0];
  int k = values[1];
  int t = values[2];
  bool qc_block_size = values[3];

  for (int i = 0; i < static_cast<int>(Algorithm::Count); i++) {
    Algorithm algo = static_cast<Algorithm>(i);
    std::cout << "Algorithm " << algorithm_to_string(algo) << std::endl;
    Result current_c_res = c_isd_log_cost(n, k, t, qc_block_size,
                                          QCAttackType::Plain, false, {algo});
    std::cout << "Plain " << std::endl;
    std::cout << result_to_string(current_c_res) << std::endl;

    double red_fac;
    red_fac =
        get_qc_red_factor_classic_log(qc_block_size, n - k, QCAttackType::MRA);
    std::cout << "Classic MRA: " << current_c_res.value - red_fac << std::endl;
    red_fac =
        get_qc_red_factor_classic_log(qc_block_size, n - k, QCAttackType::KRA1);
    std::cout << "Classic KRA1: " << current_c_res.value - red_fac << std::endl;
    red_fac =
        get_qc_red_factor_classic_log(qc_block_size, n - k, QCAttackType::KRA2);
    std::cout << "Classic KRA2: " << current_c_res.value - red_fac << std::endl;
    red_fac =
        get_qc_red_factor_classic_log(qc_block_size, n - k, QCAttackType::KRA3);
    std::cout << "Classic KRA3: " << current_c_res.value - red_fac << std::endl;

    std::cout << "**********" << std::endl;
  }

  for (int i = 0; i < static_cast<int>(QuantumAlgorithm::Count); i++) {
    QuantumAlgorithm algo = static_cast<QuantumAlgorithm>(i);

    std::cout << "Algorithm: " << quantum_algorithm_to_string(algo)
              << std::endl;
    Result current_q_res =
        q_isd_log_cost(n, k, t, qc_block_size, QCAttackType::Plain, false,
                       std::unordered_set<QuantumAlgorithm>{algo});
    std::cout << "Plain " << std::endl;
    std::cout << result_to_string(current_q_res) << std::endl;

    double red_fac;
    red_fac =
        get_qc_red_factor_quantum_log(qc_block_size, n - k, QCAttackType::MRA);
    std::cout << "Quantum MRA: " << current_q_res.value - red_fac << std::endl;
  }
  return 0;
}

int handle_json(std::string json_filename) {
  // const std::string input_isd_values = "out/isd_values.json";
  std::ifstream file(json_filename);

  // Check if the file is open
  if (!file.is_open()) {
    std::cerr << "Could not open the input file " << json_filename << std::endl;
    return 1;
  }

  // Parse the JSON content
  nlohmann::json j;
  file >> j;

  int no_values = j.size();
  std::cout << "Number of values in the JSON: " << no_values << std::endl;
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

  // Define an atomic counter for processed entries
  std::atomic<int> processed_count(0);
  std::atomic<int> error_count(0);
  std::atomic<int> skipped_count(0);

  // Iterate over the list of entries. With schedule(dynamic) loop iterations
  // are divided into chunks, and threads dynamically grab chunks as they
  // complete their previous work.
#pragma omp parallel for schedule(dynamic)
  for (const auto &entry : j) {

    if (processed_count % 1000 == 0) {
#pragma omp critical
      {
        std::cout << "\rProcessed: " << processed_count << " / " << no_values
                  << "; Skipped:" << skipped_count
                  << "; Errors: " << error_count << std::flush;
      }
    }

    uint32_t n = entry["n"];
    uint32_t r = entry["r"];
    uint32_t k = n - r;
    uint32_t t = entry["t"];

    std::string filename =
        OUT_DIR_RESULTS + fmt::format("{:06}_{:06}_{:03}.json", n, k, t);
    // Check if the generated file exists
    if (std::filesystem::exists(filename)) {
      // std::cout << "Generated file exists: " << filename << std::endl
      //           << ". Skipping.";
      continue;
      ++skipped_count;
    }
    // uint32_t qc_block_size = entry["prime"];
    uint32_t qc_block_size = r;

    nlohmann::json out_values;

    Result current_c_res;
    Result current_q_res;

    current_c_res = c_isd_log_cost(
        n, k, t, qc_block_size, QCAttackType::Plain, false,
        std::unordered_set<Algorithm>{Algorithm::Prange, Algorithm::Stern});

    current_q_res = q_isd_log_cost(
        n, k, t, qc_block_size, QCAttackType::Plain, false,
        std::unordered_set<QuantumAlgorithm>{QuantumAlgorithm::Q_Lee_Brickell});

    std::string attack_type;
    out_values["Classic"]["Plain"] = current_c_res;
    out_values["Quantum"]["Plain"] = current_q_res;

    // Post-apply reduction factors
    double red_fac =
        get_qc_red_factor_quantum_log(qc_block_size, n - k, QCAttackType::MRA);
    out_values["Quantum"]["MRA"] = current_q_res.value - red_fac;

    red_fac =
        get_qc_red_factor_classic_log(qc_block_size, n - k, QCAttackType::MRA);
    out_values["Classic"]["MRA"] = current_c_res.value - red_fac;
    red_fac =
        get_qc_red_factor_classic_log(qc_block_size, n - k, QCAttackType::KRA1);
    out_values["Classic"]["KRA1"] = current_c_res.value - red_fac;
    red_fac =
        get_qc_red_factor_classic_log(qc_block_size, n - k, QCAttackType::KRA2);
    out_values["Classic"]["KRA2"] = current_c_res.value - red_fac;
    red_fac =
        get_qc_red_factor_classic_log(qc_block_size, n - k, QCAttackType::KRA3);
    out_values["Classic"]["KRA3"] = current_c_res.value - red_fac;

    std::ofstream file(filename);
    if (file.is_open()) {
      file << std::fixed << std::setprecision(10)
           << out_values.dump(4); // Format JSON with indentation
      file.close();
      ++processed_count;
      // std::cout << "Data written to " << filename << std::endl;
    } else {
      std::cerr << "Could not open the file!" << std::endl;
      ++error_count;
    }
    if (processed_count % 1000 == 0) {
#pragma omp critical
      { std::cout << processed_count << " / " << no_values << std::endl; }
    }
  }
  return 0;
  }

int main(int argc, char *argv[]) {
  // Logger::LoggerManager::getInstance().setup_logger(
  //     "binomials", spdlog::level::info, spdlog::level::debug);
  // Logger::LoggerManager::getInstance().setup_logger(
  //     "isd_cost_estimate", spdlog::level::info, spdlog::level::debug);
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " --json [filename] | --plain [args]"
              << std::endl;
    return 1;
  }

  NTL::RR::SetPrecision(NUM_BITS_REAL_MANTISSA);
  InitBinomials();
  pi = NTL::ComputePi_RR();

  if (strcmp(argv[1], "--json") == 0) {
    std::string json_filename = argv[2];
    handle_json(json_filename);
  } else if (strcmp(argv[1], "--plain") == 0) {
    std::string plainArgs = argv[2];
    handle_plain(plainArgs);
  } else {
    std::cerr << "Unknown argument: " << argv[1] << std::endl;
    std::cerr << "Usage: " << argv[0]
              << " --json [filename]" // "| --csv [filename]"
              << std::endl;
    return 1;
  }

  return 0;
}
