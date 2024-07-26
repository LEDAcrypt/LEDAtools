add_library(ledautils binomials.cpp logging.cpp bit_error_probabilities.cpp isd_cost_estimate.cpp)

find_library(GMP_LIB gmp)
find_library(NTL_LIB ntl)
find_library(M_LIB m)
find_package(spdlog REQUIRED)
find_package(fmt REQUIRED)

message(STATUS "gmp library: ${GMP_LIBRARIES}")
message(STATUS "ntl library: ${NTL_LIBRARIES}")
message(STATUS "m library: ${M_LIBRARIES}")
message(STATUS "spdlog library: ${spdlog_LIBRARIES}")
message(STATUS "fmt library: ${fmt_LIBRARIES}")

# Specify include directories for this module
target_include_directories(ledautils PUBLIC ${PROJECT_SOURCE_DIR}/include/utils)
