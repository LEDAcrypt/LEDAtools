cmake_minimum_required(VERSION 3.10)

# Project name and version
project(LEDAtools VERSION 1.0 LANGUAGES CXX)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Global compiler flags
add_compile_options(-O3 -g3 -Wall -Wextra -Wno-sign-compare)

# Include directories
include_directories(include)

# Add the src directory
add_subdirectory(src)

# # Add the example directory
# add_subdirectory(examples)

# # Add the test directory
# add_subdirectory(test)
