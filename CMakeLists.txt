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
set(LIBS ntl gmp m )

# Include directories (if any)
# include_directories(${CMAKE_SOURCE_DIR}/include)

# Libraries
find_package(OpenMP)
find_library(gmp gmp)
find_library(ntl ntl)
find_library(m m)

find_package(spdlog QUIET)
# If spdlog is not found, add the local include directory
if (NOT spdlog_FOUND)
    message(WARNING "spdlog not found, using local include directory at $ENV{HOME}/vc/spdlog/include")
    include_directories($ENV{HOME}/vc/spdlog/include)
    link_directories($ENV{HOME}/vc/spdlog/build)
    list(APPEND LIBS libspdlog.a)
else()
    list(APPEND LIBS spdlog)
    set(LIBNAMES "spdlog")
endif()

find_package(fmt QUIET)
if (NOT fmt_FOUND)
    message(WARNING "fmt not found, using local include directory at $ENV{HOME}/vc/fmt/include")
    include_directories($ENV{HOME}/vc/fmt/include)
    link_directories($ENV{HOME}/vc/fmt/build)
    list(APPEND LIBS libfmt.a)
else()
    message(STATUS "fmt found: ${fmt_DIR}")
    list(APPEND LIBS fmt)
endif()


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
