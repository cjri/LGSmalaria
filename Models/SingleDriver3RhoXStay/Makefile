CC	      = g++
CC_FLAGS	= -g3 -O3 -Wall -I  /Users/ci3/local/include/gsl/
LD_FLAGS	= -L/Users/ci3/local/include/gsl/lib/ -lm -lgsl -lcblas
ONE		= onetimesel3rhoXs.o utilities.o

one : $(ONE)
	$(CC) $(CC_FLAGS) $(ONE) -o run_onetimesel3rhoXs $(LD_FLAGS)
onetimesel3rhoXs.o: onetimesel3rhoXs.cpp
	$(CC) $(CC_FLAGS) -c onetimesel3rhoXs.cpp
utilities.o: utilities.cpp
	$(CC) $(CC_FLAGS) -c utilities.cpp

