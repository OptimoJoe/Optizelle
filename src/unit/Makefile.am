AUTOMAKE_OPTIONS = gnu
EXTRA_PROGRAMS = run 
CLEANFILES = run 
AM_CPPFLAGS = -I.. ${BOOST_CPPFLAGS}
AM_CXXFLAGS = $(OPENMP_CXXFLAGS)
run_SOURCES = derivative_check.cpp quadratic_cone.cpp sdp_cone.cpp \
               capture_and_release.cpp linear_cone.cpp rosenbrock.cpp \
               run.cpp
run_LDADD = -L../peopt -lpeopt ${LAPACK_LIBS} ${BLAS_LIBS} \
             ${BOOST_UNIT_TEST_FRAMEWORK_LIB}
unit:
if HAVE_BOOST 
if HAVE_LAPACK
	make run 
else
	@echo "We require LAPACK in order to compile the unit tests."
endif
else
	@echo "We require Boost in order to compile the unit tests."
endif
.PHONY: unit 
