cmake_minimum_required(VERSION 3.10)

# Project name and version
project(MyProject VERSION 1.0 LANGUAGES CXX)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Compiler options
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g3 -Wall -Wextra -Wno-sign-compare")

# Define the executable targets
set(TARGETS 
    constant_weight_encodable_bits
    enumeration_complexity
    parameter_generator
    work_factor_computation
)
# Define libraries
set(LIBS ntl gmp m spdlog fmt)

# Include directories (if any)
# include_directories(${CMAKE_SOURCE_DIR}/include)

# Libraries
find_library(gmp gmp)
find_library(ntl ntl)
find_library(m m)
find_package(OpenMP REQUIRED)
find_package(spdlog REQUIRED)
find_package(fmt REQUIRED)

# Print include directories to verify
get_directory_property(dirs INCLUDE_DIRECTORIES)
message(STATUS "Include directories: ${dirs}")
message(STATUS "Libraries: ${LIBS}")

# Define the output directory
set(BIN_DIR ${CMAKE_BINARY_DIR}/bin)
file(MAKE_DIRECTORY ${BIN_DIR})

foreach(target ${TARGETS})
    add_executable(${target} ${target}.cpp)
    set_target_properties(${target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})
    target_link_libraries(${target} ${LIBS})
endforeach()

# Special case for work_factor_computation_parallel with OpenMP
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
list(APPEND LIBS OpenMP::OpenMP_CXX)
add_executable(work_factor_computation_parallel work_factor_computation_parallel.cpp)
set_target_properties(work_factor_computation_parallel PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})
target_link_libraries(work_factor_computation_parallel ${LIBS})
