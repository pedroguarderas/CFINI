GCC = g++
RM = rm
LN = ln
FLIBS = -llapack -lblas
FFLAGS = -c -g -Wall -std=c++0x -m64 -fopenmp -fPIC 
FFLAGSOPT = -O3 -fopenmp
DISLIN_INC = /usr/local/dislin/df/real64


all:
	$(MAKE) exo_1

exo_1:
	$(GCC) $(FFLAGS) exo_1.cpp 
	$(GCC) -o exo_1 exo_1.o
	
cl_exo_1:
	$(RM) exo_1.o exo_1
	
exo_2:
	$(GCC) $(FFLAGS) exo_2.cpp 
	$(GCC) -o exo_2 exo_2.o
	
cl_exo_2:
	$(RM) exo_2.o exo_2
