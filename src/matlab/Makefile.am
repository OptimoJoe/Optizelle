AUTOMAKE_OPTIONS = gnu
EXTRA_DIST = getProblemClass.cpp 
MEXFILES = getProblemClass.$(MEXEXT)
SUFFIXES = .$(MEXEXT)
pkgdata_DATA = ${MEXFILES}
.cpp.$(MEXEXT) : 
	$(MEX) $(MEXOPTS) $(MEXFLAGS) -I.. $< ../peopt/.libs/libpeopt.a \
            ${BLAS_LIBS} $(MEXLDADD)
all: ${MEXFILES}
mostlyclean-local:
	test -z '$(MEXEXT)' || rm -f *.$(MEXEXT)

