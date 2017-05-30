CC		= g++
CC_FLAGS        = -Wall -O3 -I/Users/ci3/local/include/gsl -I/opt/local/include/ -DHAVE_INLINE -fopenmp
LD_FLAGS        = -L/Users/ci3/local/include/gsl/lib/ -lm -lgsl -lcblas -fopenmp
EMISSION	= emission
JUMPDIFF	= jump-diffusion
FILTER		= jump-diffusion-filter
OBJECTS      	= $(EMISSION).o $(JUMPDIFF).o  $(FILTER).o
all:	$(OBJECTS)
	$(CC) $(CC_FLAGS) $(EMISSION).o $(JUMPDIFF).o $(FILTER).o -o ./build/$(FILTER) $(LD_FLAGS)
	rm -f ./*.o
$(FILTER).o: $(FILTER).cpp
	$(CC) $(CC_FLAGS) -c $(FILTER).cpp
$(EMISSION).o: $(EMISSION).cpp $(EMISSION).h
	$(CC) $(CC_FLAGS) -c $(EMISSION).cpp 
$(JUMPDIFF).o: $(JUMPDIFF).cpp $(JUMPDIFF).h
	$(CC) $(CC_FLAGS) -c $(JUMPDIFF).cpp  
clean:
	rm -f ./*.o

