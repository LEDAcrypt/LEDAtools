CXXFLAGS=-O3 -g3 -std=c++11 -Wall -Wextra -Wno-sign-compare
OMPFLAGS=-fopenmp
LDLIBS= -lntl -lgmp -lm
BIN_DIR=bin
TARGETS=constant_weight_encodable_bits enumeration_complexity parameter_generator work_factor_computation work_factor_computation_parallel
BIN_TARGETS=$(addprefix $(BIN_DIR)/, $(TARGETS))

all: $(BIN_DIR) $(BIN_TARGETS)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(BIN_DIR)/work_factor_computation_parallel: work_factor_computation_parallel.cpp
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) $< -o $@ $(LDLIBS)

$(BIN_DIR)/%: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDLIBS)

clean:
	rm -f $(BIN_DIR)/*
	rmdir $(BIN_DIR)

