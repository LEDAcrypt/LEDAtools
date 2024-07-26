
find_package(OpenMP REQUIRED)
message(STATUS "OpenMP library: ${OpenMP_CXX_LIBRARIES}")
find_library(GMP_LIB gmp)
find_library(NTL_LIB ntl)
find_library(M_LIB m)
find_package(spdlog REQUIRED)
find_package(fmt REQUIRED)

# Define libraries
set(LIBS ledautils ${NTL_LIB} ${GMP_LIB} ${M_LIB} spdlog::spdlog fmt::fmt)

message(STATUS "gmp library: ${GMP_LIBRARIES}")
message(STATUS "ntl library: ${NTL_LIBRARIES}")
message(STATUS "m library: ${M_LIBRARIES}")
message(STATUS "spdlog library: ${spdlog_LIBRARIES}")
message(STATUS "fmt library: ${fmt_LIBRARIES}")

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

install(TARGETS constant_weight_encodable_bits DESTINATION bin)
install(TARGETS enumeration_complexity DESTINATION bin)
install(TARGETS parameter_generator DESTINATION bin)
install(TARGETS work_factor_computation DESTINATION bin)
install(TARGETS work_factor_computation_parallel DESTINATION bin)
