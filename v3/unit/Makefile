include ../make.inc

# List the sources for compilation
OBJS = capture_and_release.o derivative_check.o rosenbrock.o linear_cone.o \
    quadratic_cone.o sdp_cone.o unit.o

CFLAGS := ${CFLAGS}
INCLUDES := ${INCLUDES} -I..
LDFLAGS := -L.. -lpeopt ${LDFLAGS}
      
all: ${OBJS}
	${CPP} ${OBJS} ${LDFLAGS} -o unit

%.o : %.cpp
	${CPP} ${INCLUDES} ${CFLAGS} -c $<

clean:
	rm ${OBJS} unit 
