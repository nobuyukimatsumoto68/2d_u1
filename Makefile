CXX		= g++
CXXFLAGS	= -O3 -Wall -MMD -MP -std=c++17 # -fopenmp
PROG		= a.out
SRC		= main.c

OBJ		= $(SRC:%.c=%.o)
DEPS		= $(SRC:%.c=%.d)
-include $(DEPS)

$(OBJ): $(SRC)
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) $< -o $@

# .PHONY: clean
# clean:;	rm -f *.out
