include make.inc

# List the sources for compilation
OBJS = linalg.o

all: $(OBJS)
	ar cr libpeopt.a ${OBJS}

%.o : %.cpp
	${CPP} ${INCLUDES} ${CFLAGS} -c $<

clean:
	rm ${OBJS} libpeopt.a 
