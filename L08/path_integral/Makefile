LIBS:=`root-config --libs`
INCS:=`root-config --cflags`
exe = qmc1d_variational qmc1d_uniform

all : $(exe)

qmc1d_uniform: qmc1d_uniform.o
	g++ -O3 -Wall -o $@ $^ ${LIBS}
	
qmc1d_variational: qmc1d_variational.o
	g++ -O3 -Wall -o $@ $^ ${LIBS}

%.o : %.cpp
	g++ -Wall -c $< ${INCS}

clean:
	rm -rf *.o $(exe)
