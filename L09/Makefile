CC = g++
CFLAGS = -Wall -g -O3 --std=c++11

head =  random.h
obj =  random.o
exe  =  exe-09
#do not set more than one executable!

.PHONY: clean

all: $(exe)

$(exe) : $(obj) $(exe).o
	$(CC) $^  -o $(@)
	
%.o: %.cpp $(head)
	$(CC) $(CFLAGS) -c $< -o $@
	
clean :
	rm *.o $(exe)
