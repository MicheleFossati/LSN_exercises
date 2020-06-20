CC = g++
CFLAGS = -Wall -g -O3 --std=c++11

head =  random.h
obj =  random.o
exe  =  exe03

.PHONY: clean

all: $(exe)

$(exe) : $(obj) $(exe).o
	$(CC) $^  -o $(@)
	
%.o: %.cpp $(head)
	$(CC) $(CFLAGS) -c $< -o $@
	
clean :
	rm *.o $(exe)

