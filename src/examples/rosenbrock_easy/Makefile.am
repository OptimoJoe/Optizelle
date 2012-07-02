AUTOMAKE_OPTIONS = gnu
EXTRA_PROGRAMS = run 
CLEANFILES = run 
AM_CPPFLAGS = -I../.. ${JSONCPP_CPPFLAGS}
AM_CXXFLAGS = $(OPENMP_CXXFLAGS)
run_SOURCES = run.cpp
run_LDADD = -L../../peopt -lpeopt ${LAPACK_LIBS} ${BLAS_LIBS} \
             ${JSONCPP_LDFLAGS} ${JSONCPP_LIBS}
example:
if HAVE_JSONCPP
if HAVE_LAPACK
	make run 
else
	@echo "We require LAPACK in order to compile this example."
endif
else
	@echo "We require jsoncpp in order to compile this example."
endif
.PHONY: example 
