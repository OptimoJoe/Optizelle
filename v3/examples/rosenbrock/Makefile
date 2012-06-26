include ../../make.inc

OBJS = run.o
INCLUDES := ${INCLUDES} -I../..

all: $(OBJS)
	${CPP} ${LDFLAGS} ${OBJS} -o run

%.o : %.cpp
	${CPP} ${INCLUDES} ${CFLAGS} -c $<

clean:
	rm ${OBJS} run 
