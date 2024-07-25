# Define the executable targets
set(TARGETS 
    constant_weight_encodable_bits
    enumeration_complexity
    parameter_generator
    work_factor_computation
    work_factor_computation_parallel
)

# Libraries
# Find libraries
find_library(GMP_LIB gmp)
find_library(NTL_LIB ntl)
find_library(M_LIB m)
find_package(OpenMP REQUIRED)
find_package(spdlog REQUIRED)
find_package(fmt REQUIRED)

message(STATUS "gmp library: ${GMP_LIBRARIES}")
message(STATUS "ntl library: ${NTL_LIBRARIES}")
message(STATUS "m library: ${M_LIBRARIES}")
message(STATUS "OpenMP library: ${OpenMP_CXX_LIBRARIES}")
message(STATUS "spdlog library: ${spdlog_LIBRARIES}")
message(STATUS "fmt library: ${fmt_LIBRARIES}")

# Define libraries
set(LIBS ${NTL_LIB} ${GMP_LIB} ${M_LIB} spdlog::spdlog fmt::fmt)

foreach(target ${TARGETS})
    add_executable(${target} ${target}.cpp)
    target_link_libraries(${target} PRIVATE ${LIBS})
endforeach()

# Special case for work_factor_computation_parallel with OpenMP
target_link_libraries(work_factor_computation_parallel PRIVATE OpenMP::OpenMP_CXX)

# message(STATUS "m library: ${M_LIBRARIES}")
