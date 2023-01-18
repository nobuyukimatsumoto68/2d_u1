CXX		= g++
CXXFLAGS	= -O3 -Wall -MMD -MP -std=c++17 -fopenmp
# CXXFLAGS	= -O1 -Wall -MMD -MP -std=c++17 -pg
PROG		= a.out
SRC		= main.c

OBJ		= $(SRC:%.c=%.o)
DEPS		= $(SRC:%.c=%.d)
-include $(DEPS)

$(OBJ): $(SRC)
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) $< -o $@

.PHONY: run
run:; time ./main.o 2>&1 | tee log

# .PHONY: clean
# clean:;	rm *.out
