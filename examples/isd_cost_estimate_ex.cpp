#include "isd_cost_estimate.hpp"
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
  uint32_t n0;
  uint32_t prime;
  uint32_t v;
  uint32_t t;
  std::vector<Cost> costs;
};

void displayValues(const std::vector<Value> &values) {
  // Optional: Display the values
  for (const auto &value : values) {
    std::cout << "n0: " << value.n0 << ", prime: " << value.prime << "\n";
    for (const auto &cost : value.costs) {
      std::cout << "  Algorithm: " << cost.algorithm << ", Type: " << cost.type
                << ", Quantum: " << (cost.is_quantum ? "Yes" : "No")
                << ", Time Complexity: " << cost.time_complexity
                << ", Space Complexity: " << cost.space_complexity << "\n";
    }
  }
}

int main() {
  std::cout << "Hello world\n";
  // Expected values taken from LEDA specs, Table 4.1
  std::vector<Value> values;
  Value val = {2,     
               23371, 
               71,    
               130,   
               {
                   {
                       "Prange", 
                       "CFP1",   
                       false,    
                       144.2,    
                       0.0       
                   },
                   {
                       "Prange",
                       "CFP1",  
                       false,   
                       144.2,   
                       0.0      
                   },

               }};

  // Value value1 = {
  //         10, // n0
  //         7,  // prime
  //         5,  // v
  //         20, // t
  //         {{"Algorithm1", "CFP1", true, 0.5, 1.0},
  //          {"Algorithm2", "SDP", false, 1.0, 2.0}} // costs
  //     };
  // values.push_back(value1);
  // ,
  //   {
  //       15, // n0
  //       11, // prime
  //       6,  // v
  //       30, // t
  //       {{"Algorithm3", "CFP2", true, 0.7, 1.5},
  //        {"Algorithm4", "CFP3", false, 0.8, 1.8}} // costs
  //   }};

  // Call the function to display the values
  // displayValues(values);

  // return 0;
}
