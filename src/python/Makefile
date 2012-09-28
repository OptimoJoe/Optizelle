INCLUDE = -I/home/josyoun/usr/include/python2.7 -I../peopt
LDFLAGS = -L/home/josyoun/usr/lib/python2.7/config -lpython2.7 -lpthread -ldl -lutil
CFLAGS = -fPIC -g
.cpp.o:
	g++ -c ${CFLAGS} ${INCLUDE} $< 
all: peopt.o 
	g++ -shared -Wl,-soname,libpeopt.so -o libpeopt.so peopt.o 
clean:
	rm peopt.o libpeopt.so peopt.pyc
