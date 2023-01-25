SRC		= main.cpp
OBJ		= $(SRC:%.cpp=%.o)
DIR_SRC		= src
IS_FTHMC	= 0

export SRC OBJ IS_FTHMC

all:
	@make -C $(DIR_SRC)
	@cp $(DIR_SRC)/$(OBJ) .

.PHONY: run
run:; time ./$(OBJ) 2>&1 | tee log
