CC = g++
MPICC = mpic++
CFLAGS = -Wall -g -O3 --std=c++11

head =  random.h
obj =  random.o
exe  =  exe10-2

.PHONY: clean 

all: $(exe)

$(exe) : $(obj) $(exe).o
	$(MPICC) $^  -o $(@)
	
$(exe).o : $(exe).cpp
	$(MPICC) $(CFLAGS) -c $< -o $@

%.o: %.cpp $(head)
	$(CC) $(CFLAGS) -c $< -o $@

exe : 
	mpiexec -np 4 ./$(exe)	
clean :
	rm *.o $(exe)

load :
	module load mpi/mpich-3.2-x86_64

