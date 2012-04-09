# List the sources for compilation
OBJS = test01.o

# Set the optimizations
CFLAGS= -Wall -g #-O2

all: $(OBJS)
	g++ test01.o -o test01

%.o : %.cpp
	g++ $(CFLAGS) -c $<

clean:
	rm *.o test01
