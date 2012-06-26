include ../../make.inc

OBJS = run.o
INCLUDES := ${INCLUDES} -I../..
LDFLAGS := -L../.. -lpeopt ${LDFLAGS}

all: $(OBJS)
	${CPP} ${OBJS} ${LDFLAGS} -o run

%.o : %.cpp
	${CPP} ${INCLUDES} ${CFLAGS} -c $<

clean:
	rm ${OBJS} run 
