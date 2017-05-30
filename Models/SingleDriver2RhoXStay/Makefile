CC	      = g++
CC_FLAGS	= -g3 -O3 -Wall -I  /Users/ci3/local/include/gsl/
LD_FLAGS	= -L/Users/ci3/local/include/gsl/lib/ -lm -lgsl -lcblas
ONE		= onetimesel2rhoXs.o utilities.o

one : $(ONE)
	$(CC) $(CC_FLAGS) $(ONE) -o run_onetimesel2rhoXs $(LD_FLAGS)
onetimesel2rhoXs.o: onetimesel2rhoXs.cpp
	$(CC) $(CC_FLAGS) -c onetimesel2rhoXs.cpp
utilities.o: utilities.cpp
	$(CC) $(CC_FLAGS) -c utilities.cpp

