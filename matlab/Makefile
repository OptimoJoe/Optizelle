INCLUDES=-I/home/josyoun/usr/include -I/opt/matlab/extern/include/
CFLAGS=-g -fPIC
LDFLAGS=-L/home/josyoun/usr/lib -llapack -lf77blas -lcblas -latlas
OBJECTS=test01.o simple_matching.o

.cpp.o:
	g++ -c ${CFLAGS} ${INCLUDES} $<

all: ${OBJECTS} 
	g++ ${OBJECTS} ${LDFLAGS} -o test01

clean:
	rm *.o *.mexa64 test01
