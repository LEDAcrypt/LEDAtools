# Define the executable targets
set(TARGETS 
    constant_weight_encodable_bits
    enumeration_complexity
    parameter_generator
    work_factor_computation
    work_factor_computation_parallel
)

# Libraries
find_library(gmp gmp)
find_library(ntl ntl)
find_library(m m)
# Packages
message(STATUS "m library: ${M_LIBRARIES}")
find_package(OpenMP REQUIRED)
message(STATUS "OpenMP library: ${OpenMP_CXX_LIBRARIES}")
find_package(spdlog REQUIRED)
message(STATUS "spdlog library: ${spdlog_LIBRARIES}")
find_package(fmt REQUIRED)
message(STATUS "fmt library: ${fmt_LIBRARIES}")

# Define libraries
set(LIBS ntl gmp m spdlog fmt)

foreach(target ${TARGETS})
    add_executable(${target} ${target}.cpp)
    target_link_libraries(${target} PRIVATE ${LIBS})
endforeach()

# Special case for work_factor_computation_parallel with OpenMP
target_link_libraries(work_factor_computation_parallel PRIVATE OpenMP::OpenMP_CXX)
