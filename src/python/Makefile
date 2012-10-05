INCLUDE = -I/home/josyoun/usr/include/python2.7 -I.. -I/home/josyoun/usr/include
LDFLAGS = ../peopt/.libs/libpeopt.a -L/home/josyoun/usr/lib -ljson
CFLAGS = -fPIC -g
.cpp.o:
	g++ -c ${CFLAGS} ${INCLUDE} $< 
all: peopt.o 
	g++ -shared -o libpeopt.so peopt.o ${LDFLAGS}
clean:
	rm peopt.o libpeopt.so peopt.pyc
