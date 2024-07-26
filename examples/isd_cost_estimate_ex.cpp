#include <isd_cost_estimate.hpp>
#include <cstdint> // for uint32_t
#include <vector>
#include <string>
#include <ostream>
#include <iostream>

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
  std::vector<Cost> costs;
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

  std::cout << "Costs:\n";
  for (const auto &cost : value.costs) {
    displayCost(cost);
    std::cout << "-----\n";
  }
}

int main() {
  std::cout << "Hello world\n";
  // Expected values taken from LEDA specs, Table 4.1
  std::vector<Value> values;
  Value val = {24646,
               12323,
               142,
               12323,
               true,
               {
                   {"Prange", "", false, 171.3, 0.0},
                   {"Lee-Brickell", "", false, 158.4, 0.0},
                   {"Leon", "", false, 154.4, 0.0},
                   {"Stern", "", false, 147.4, 0.0},
                   {"Fin-Send", "", false, 147.4, 0.0},
               }};

  displayValue(val);
}
