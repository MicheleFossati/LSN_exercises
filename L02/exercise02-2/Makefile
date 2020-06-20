CC = g++
CFLAGS = -Wall -g -O3 --std=c++11

head =  random.h walker.h
obj =  random.o walker.o
exe  =  exe02-2

.PHONY: clean

all: $(exe)

$(exe) : $(obj) $(exe).o
	$(CC) $^  -o $(@)
	
%.o: %.cpp $(head)
	$(CC) $(CFLAGS) -c $< -o $@
	
clean :
	rm *.o $(exe)

