cmake_minimum_required(VERSION 3.10)

# Project name and version
project(LEDAtools VERSION 1.0 LANGUAGES CXX)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Global compiler flags
add_compile_options(-O3 -g3 -Wall -Wextra -Wno-sign-compare)

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

add_subdirectory(src)
add_subdirectory(examples)
# add_subdirectory(test)

# Installation directories
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
