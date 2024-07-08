CXXFLAGS=-O3 -g3 -std=c++11 -Wall -Wextra -Wno-sign-compare
LDLIBS= -lntl -lgmp -lm
TARGETS=constant_weight_encodable_bits enumeration_complexity parameter_generator work_factor_computation

all: $(TARGETS)

clean:
	rm -f $(TARGETS)
