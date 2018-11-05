##
## on mahuika: sbatch run.sl
##   check using: squeue, sacct
##
###########################################################################################
## C++
CPP_DEFS = -std=c++11
CPP = mpicxx

###########################################################################################
## NOTE: uncomment block for ONE of debug, profile or optimised

## debug
#CPP_FLAGS =  $(CPP_DEFS) $(CPP_INCLUDES) -m64 -Og -Wall -Wextra

## profile
#CPP_FLAGS =  $(CPP_DEFS) $(CPP_INCLUDES) -m64 -O3 -pg 

## optimised
CPP_FLAGS =  $(CPP_DEFS) $(CPP_INCLUDES) -m64 -O3 -march=native -mtune=native 

###########################################################################################
## program

## Mahuika
CPP_INCLUDES = -I/usr/include -I./src
## Workstation
#CPP_INCLUDES = -I/usr/include -I/usr/include/eigen3

EXE = ../run/parotid_sim_5
OBJS = main.o utils.o cAcinus.o cCell_x.o cVCLSolver.o cCellMesh.o
SRC = ./src

all: $(EXE)

main.o: $(SRC)/main.cpp $(SRC)/global_defs.hpp
	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

utils.o: $(SRC)/utils.cpp $(SRC)/utils.hpp
	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

cAcinus.o: $(SRC)/cAcinus.cpp $(SRC)/cAcinus.hpp
	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

cCell_x.o: $(SRC)/cCell_x.cpp $(SRC)/cCell_x.hpp
	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

cVCLSolver.o: $(SRC)/cVCLSolver.cpp $(SRC)/cVCLSolver.hpp
	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

cCellMesh.o: $(SRC)/cCellMesh.cpp $(SRC)/cCellMesh.hpp
	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

$(EXE): $(OBJS)
	$(CPP) $(CPP_FLAGS) -o $@ $+

clean:
	rm *.o


