CC              = g++
CC_FLAGS        = -g3 -O3 -Wall -I  /Users/ci3/local/include/gsl/
LD_FLAGS        = -L/Users/ci3/local/include/gsl/lib/ -lm -lgsl -lcblas
ONE		= onetimesel.o utilities.o

one : $(ONE)
	$(CC) $(CC_FLAGS) $(ONE) -o run_onetimesel $(LD_FLAGS)
onetimesel.o: onetimesel.cpp
	$(CC) $(CC_FLAGS) -c onetimesel.cpp
utilities.o: utilities.cpp
	$(CC) $(CC_FLAGS) -c utilities.cpp

show : $(SHOW)
	$(CC) $(CC_FLAGS) $(ONE) -o showSD $(LD_FLAGS)
showSD.o: showSD.cpp
	$(CC) $(CC_FLAGS) -c showSD.cpp

two : $(TWO)
	$(CC) $(CC_FLAGS) $(TWO) -o run_onetimesel_range $(LD_FLAGS)
onetimesel_range.o: onetimesel_range.cpp
	$(CC) $(CC_FLAGS) -c onetimesel_range.cpp
utilities_range.o: utilities_range.cpp
	$(CC) $(CC_FLAGS) -c utilities_range.cpp

