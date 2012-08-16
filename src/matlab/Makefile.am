AUTOMAKE_OPTIONS = gnu
EXTRA_DIST = peopt.cpp 
MEXFILES = peopt.$(MEXEXT)
SUFFIXES = .$(MEXEXT)
pkgdata_DATA = ${MEXFILES}
.cpp.$(MEXEXT) : 
	$(MEX) $(MEXOPTS) $(MEXFLAGS) -I.. ${JSONCPP_CPPFLAGS} $< \
            ../peopt/.libs/libpeopt.a \
            ${BLAS_LIBS} ${JSONCPP_LDFLAGS} ${JSONCPP_LIBS} $(MEXLDADD)
all: ${MEXFILES}
mostlyclean-local:
	test -z '$(MEXEXT)' || rm -f *.$(MEXEXT)
