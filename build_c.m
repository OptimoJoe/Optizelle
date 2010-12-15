% Builds the c interface to the trust region code
function build_c()

% Make sure all of the object files are up to date
system('make');

% Then, build the interfaces
mex -g -I/home/josyoun/usr/include -L/home/josyoun/usr/lib -llapack -lf77blas -lcblas -latlas f.cpp simple_matching.o
mex -g -I/home/josyoun/usr/include -L/home/josyoun/usr/lib -llapack -lf77blas -lcblas -latlas fp.cpp simple_matching.o
mex -g -I/home/josyoun/usr/include -L/home/josyoun/usr/lib -llapack -lf77blas -lcblas -latlas fps.cpp simple_matching.o
mex -g -I/home/josyoun/usr/include -L/home/josyoun/usr/lib -llapack -lf77blas -lcblas -latlas fpps.cpp simple_matching.o
mex -g -I/home/josyoun/usr/include -L/home/josyoun/usr/lib -llapack -lf77blas -lcblas -latlas h.cpp simple_matching.o
mex -g -I/home/josyoun/usr/include -L/home/josyoun/usr/lib -llapack -lf77blas -lcblas -latlas hp.cpp simple_matching.o
mex -g -I/home/josyoun/usr/include -L/home/josyoun/usr/lib -llapack -lf77blas -lcblas -latlas hps.cpp simple_matching.o
mex -g -I/home/josyoun/usr/include -L/home/josyoun/usr/lib -llapack -lf77blas -lcblas -latlas hpps.cpp simple_matching.o
mex -g -I/home/josyoun/usr/include -L/home/josyoun/usr/lib -llapack -lf77blas -lcblas -latlas getGradient.cpp simple_matching.o
mex -g -I/home/josyoun/usr/include -L/home/josyoun/usr/lib -llapack -lf77blas -lcblas -latlas GaussNewton.cpp simple_matching.o
mex -g -I/home/josyoun/usr/include -L/home/josyoun/usr/lib -llapack -lf77blas -lcblas -latlas Newton.cpp simple_matching.o
mex -g -I/home/josyoun/usr/include -L/home/josyoun/usr/lib -llapack -lf77blas -lcblas -latlas BFGSinv.cpp simple_matching.o
mex -g -I/home/josyoun/usr/include -L/home/josyoun/usr/lib -llapack -lf77blas -lcblas -latlas BFGS.cpp simple_matching.o
mex -g -I/home/josyoun/usr/include -L/home/josyoun/usr/lib -llapack -lf77blas -lcblas -latlas SR1.cpp simple_matching.o
mex -g -I/home/josyoun/usr/include -L/home/josyoun/usr/lib -llapack -lf77blas -lcblas -latlas getStep.cpp simple_matching.o
mex -g -I/home/josyoun/usr/include -L/home/josyoun/usr/lib -llapack -lf77blas -lcblas -latlas pe_test.cpp simple_matching.o
