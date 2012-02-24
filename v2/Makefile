# Set the compiler
MEX=/opt/matlab/bin/mex

# List the sources for compilation
MOBJS = pelab.mex64a gradientCheck.mex64a hessianCheck.mex64a getObjective.mex64a getGradient.mex64a fpCheck.mex64a getNewton.mex64a getGaussNewton.mex64a

# Set the optimizations
CFLAGS= -g #-O2

all: $(MOBJS)

%.mex64a : %.cpp
	$(MEX) CXXOPTIMFLAGS=$(CFLAGS) $<

clean:
	rm *.mexa64
