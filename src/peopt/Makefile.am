AUTOMAKE_OPTIONS = gnu
if HAVE_LAPACK
    lib_LTLIBRARIES = libpeopt.la
else
    lib_LTLIBRARIES =
endif
AM_CXXFLAGS = $(OPENMP_CXXFLAGS)
libpeopt_la_SOURCES = linalg.cpp json.h linalg.h peopt.h vspaces.h
if HAVE_LAPACK
    pkginclude_HEADERS = json.h linalg.h peopt.h vspaces.h
else
    pkginclude_HEADERS = json.h peopt.h
endif
