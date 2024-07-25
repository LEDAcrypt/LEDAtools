# Libraries
# Find libraries
find_library(GMP_LIB gmp)
find_library(NTL_LIB ntl)
find_library(M_LIB m)
find_package(OpenMP REQUIRED)
find_package(spdlog REQUIRED)
find_package(fmt REQUIRED)

add_library(leda_tools binomials.cpp logging.cpp bit_error_probabilities.cpp isd_cost_estimate.cpp)

message(STATUS "gmp library: ${GMP_LIBRARIES}")
message(STATUS "ntl library: ${NTL_LIBRARIES}")
message(STATUS "m library: ${M_LIBRARIES}")
message(STATUS "OpenMP library: ${OpenMP_CXX_LIBRARIES}")
message(STATUS "spdlog library: ${spdlog_LIBRARIES}")
message(STATUS "fmt library: ${fmt_LIBRARIES}")

# Define libraries
set(LIBS leda_tools ${NTL_LIB} ${GMP_LIB} ${M_LIB} spdlog::spdlog fmt::fmt)

set(target constant_weight_encodable_bits)
add_executable(${target} ${target}.cpp)
target_link_libraries(${target} PRIVATE ${LIBS})

set(target enumeration_complexity)
add_executable(${target} ${target}.cpp)
target_link_libraries(${target} PRIVATE ${LIBS})

set(target parameter_generator)
add_executable(${target} ${target}.cpp)
target_link_libraries(${target} PRIVATE ${LIBS})

set(target work_factor_computation)
add_executable(${target} ${target}.cpp)
target_link_libraries(${target} PRIVATE ${LIBS})

set(target work_factor_computation_parallel)
add_executable(${target} ${target}.cpp)
target_link_libraries(${target} PRIVATE ${LIBS} OpenMP::OpenMP_CXX)
