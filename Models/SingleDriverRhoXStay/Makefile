CC	      = g++
CC_FLAGS	= -g3 -O3 -Wall -I  /Users/ci3/local/include/gsl/
LD_FLAGS	= -L/Users/ci3/local/include/gsl/lib/ -lm -lgsl -lcblas
ONE		= onetimeselrhoXs.o utilities.o

one : $(ONE)
	$(CC) $(CC_FLAGS) $(ONE) -o run_onetimeselrhoX_stay $(LD_FLAGS)
onetimeselrhoXs.o: onetimeselrhoXs.cpp
	$(CC) $(CC_FLAGS) -c onetimeselrhoXs.cpp
utilities.o: utilities.cpp
	$(CC) $(CC_FLAGS) -c utilities.cpp

