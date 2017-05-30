CC              = g++
CC_FLAGS        = -g3 -O3 -Wall -I  /Users/ci3/local/include/gsl/
LD_FLAGS        = -L/Users/ci3/local/include/gsl/lib/ -lm -lgsl -lcblas
ONE		= onefreq.o utilities.o

one : $(ONE)
	$(CC) $(CC_FLAGS) $(ONE) -o run_onefreq $(LD_FLAGS)
onefreq.o: onefreq.cpp
	$(CC) $(CC_FLAGS) -c onefreq.cpp
utilities.o: utilities.cpp
	$(CC) $(CC_FLAGS) -c utilities.cpp
