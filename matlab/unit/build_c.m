% Builds the c interface to the trust region code
function build_c()

% Then, build the interfaces
mex -g -I/home/josyoun/usr/include -I/home/josyoun/peopt -L/home/josyoun/usr/lib -llapack -lf77blas -lcblas -latlas parest.cpp ../simple_matching.o
