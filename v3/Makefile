# List the sources for compilation
OBJS = test01.o test02.o test03.o linalg.o test04.o test05.o test06.o test07.o\
       test08.o

# Included files
INCLUDES = -I /home/josyoun/usr/include

# Libraries
LDFLAGS = -L /home/josyoun/usr/lib -llapack -lcblas -lf77blas -latlas \
    -fopenmp -lgfortran -ljson

# Set the optimizations
CFLAGS= -Wall -g -pedantic ${INCLUDES} -fopenmp #-O2

# Compiler
CPP = g++44

all: $(OBJS)
	${CPP} test01.o -o test01
	${CPP} test02.o -o test02
	${CPP} test03.o linalg.o ${LDFLAGS} -o test03
	${CPP} test04.o linalg.o ${LDFLAGS} -o test04
	${CPP} test05.o linalg.o ${LDFLAGS} -o test05
	${CPP} test06.o linalg.o ${LDFLAGS} -o test06
	${CPP} test07.o linalg.o ${LDFLAGS} -o test07
	${CPP} test08.o linalg.o ${LDFLAGS} -o test08

%.o : %.cpp
	${CPP} $(CFLAGS) -c $<

clean:
	rm *.o test01 test02 test03 test04 test05 test06 test07 test08
