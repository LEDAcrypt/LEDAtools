find_library(ledautils ledautils)
message(STATUS "ledautils library: ${ledautils_LIBRARIES}")

find_package(spdlog REQUIRED)
find_package(fmt REQUIRED)

add_executable(isd_cost_estimate_ex isd_cost_estimate_ex.cpp)
set(LIBS ledautils ${NTL_LIB} ${GMP_LIB} ${M_LIB} spdlog::spdlog fmt::fmt)
target_link_libraries(isd_cost_estimate_ex PRIVATE ${LIBS})

# # Optionally, specify include directories for the examples
# target_include_directories(isd_cost_estimate_ex PRIVATE ${CMAKE_SOURCE_DIR}/include)

install(TARGETS isd_cost_estimate_ex DESTINATION bin)
