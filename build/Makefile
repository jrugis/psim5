###########################################################################################
##
## target platform: mahuika
##
## assumes:
##   ml gimkl/2017a
##   ml Eigen/3.3.4
##   ml SUNDIALS/4.0.0-gimkl-2017a
##
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
CPP_FLAGS =  $(CPP_DEFS) $(CPP_INCLUDES) -m64 -O3 -march=native -mtune=native -Wall -g

###########################################################################################
## program

## Mahuika
#CPP_INCLUDES = -I/usr/include
CPP_LIBS = -lsundials_cvode -lsundials_nvecserial

EXE = ../run/psim5
OBJS = main.o utils.o cAcinus.o cCell_calcium.o cCellMesh.o cLumen.o LSODA.o cLSODA.o cCVode.o
SRC = ./src

all: $(EXE)

LSODA.o: $(SRC)/libsoda/LSODA.cpp $(SRC)/libsoda/LSODA.h $(SRC)/libsoda/helper.h
	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

cLSODA.o: $(SRC)/cLSODA.cpp $(SRC)/cLSODA.hpp $(SRC)/libsoda/LSODA.h $(SRC)/cLumen.hpp
	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

main.o: $(SRC)/main.cpp $(SRC)/global_defs.hpp $(SRC)/cAcinus.hpp $(SRC)/cCell_calcium.hpp
	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

utils.o: $(SRC)/utils.cpp $(SRC)/utils.hpp $(SRC)/global_defs.hpp
	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

cAcinus.o: $(SRC)/cAcinus.cpp $(SRC)/cAcinus.hpp $(SRC)/utils.hpp $(SRC)/global_defs.hpp $(SRC)/cLumen.hpp
	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

cLumen.o: $(SRC)/cLumen.cpp $(SRC)/cLumen.hpp $(SRC)/utils.hpp $(SRC)/global_defs.hpp $(SRC)/cCVode.hpp $(SRC)/cLSODA.hpp
	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

cCVode.o: $(SRC)/cCVode.cpp $(SRC)/cCVode.hpp $(SRC)/global_defs.hpp $(SRC)/utils.hpp $(SRC)/cLumen.hpp
	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

cCell_calcium.o: $(SRC)/cCell_calcium.cpp $(SRC)/cCell_calcium.hpp $(SRC)/cCellMesh.hpp $(SRC)/utils.hpp $(SRC)/global_defs.hpp
	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

cCellMesh.o: $(SRC)/cCellMesh.cpp $(SRC)/cCellMesh.hpp $(SRC)/cCell_calcium.hpp $(SRC)/utils.hpp $(SRC)/global_defs.hpp
	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

$(EXE): $(OBJS)
	$(CPP) $(CPP_FLAGS) $(CPP_LIBS) -o $@ $+

clean:
	rm *.o
