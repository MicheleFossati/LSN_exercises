CC = g++
CFLAGS = -Wall -g -O3 --std=c++11

head =  random.h wavefunction.h
obj =  random.o
exe  =  wavefunction

.PHONY: clean

all: $(exe)

$(exe) : $(obj) $(exe).o
	$(CC) $^  -o $(@)
	
%.o: %.cpp $(head)
	$(CC) $(CFLAGS) -c $< -o $@
	
clean :
	rm -rf *.o $(exe)
