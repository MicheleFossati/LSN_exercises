CC = g++
CFLAGS = -Wall -O3 --std=c++11


headers = random.h WavefMRT.h
objects = random.o 
exe  =  WavefMRT


all: $(exe)

WavefMRT : $(objects)  WavefMRT.o
	$(CC) WavefMRT.o  $(objects) -o $(@F)
	
	
$(objects) : $(headers)

WavefMRT.o : WavefMRT.cpp $(headers)
	$(CC) $(CFLAGS) -c WavefMRT.cpp -o WavefMRT.o
	

	
clean :
	rm *.o WavefMRT
