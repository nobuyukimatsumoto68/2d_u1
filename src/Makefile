CXX		= g++
CXXFLAGS	= -O3 -Wall -MMD -MP -std=c++17 -DIS_FTHMC=$(IS_FTHMC) # -fopenmp

DEPS		= $(SRC:%.cpp=%.d)
-include $(DEPS)

$(OBJ): $(SRC)
	$(CXX) $(CXXFLAGS) $< -o $@
