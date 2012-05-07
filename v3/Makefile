# List the sources for compilation
OBJS = test01.o test02.o

# Set the optimizations
CFLAGS= -Wall -g -pedantic #-O2

all: $(OBJS)
	g++ test01.o -o test01
	g++ test02.o -o test02

%.o : %.cpp
	g++ $(CFLAGS) -c $<

clean:
	rm *.o test01 test02
