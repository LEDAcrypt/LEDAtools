add_executable(isd_cost_estimate_ex isd_cost_estimate_ex.cpp)
target_link_libraries(isd_cost_estimate_ex PRIVATE leda_tools)
# Optionally, specify include directories for the examples
target_include_directories(isd_cost_estimate_ex PRIVATE ${CMAKE_SOURCE_DIR}/include)
