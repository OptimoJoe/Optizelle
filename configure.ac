# Initialization
AC_PREREQ(2.59)
AC_INIT(peopt, 0.3.0, josyoun@sandia.gov)
AC_CONFIG_AUX_DIR(build-aux)
AM_INIT_AUTOMAKE(peopt, 0.3.0, josyoun@sandia.gov)

# Checks for programs.
AC_PROG_CXX       # Needed for the primary program
AC_PROG_F77       # Needed for testing BLAS, LAPACK, and name mangling
AC_PROG_LIBTOOL   # Needed for libtool

# Checks for libraries.

# Checks for standard header files.
AC_LANG_PUSH(C++)
AC_CHECK_HEADERS(list,,
    [AC_MSG_ERROR([The list header is required for peopt])])
AC_CHECK_HEADERS(vector,,
    [AC_MSG_ERROR([The vector header is required for peopt])])
AC_CHECK_HEADERS(limits,,
    [AC_MSG_ERROR([The limits header is required for peopt])])
AC_CHECK_HEADERS(cmath,,
    [AC_MSG_ERROR([The cmath header is required for peopt])])
AC_CHECK_HEADERS(sstream,,
    [AC_MSG_ERROR([The sstream header is required for peopt])])
AC_CHECK_HEADERS(iostream,,
    [AC_MSG_ERROR([The iostream header is required for peopt])])
AC_CHECK_HEADERS(iomanip,,
    [AC_MSG_ERROR([The iomanip header is required for peopt])])
AC_CHECK_HEADERS(memory,,
    [AC_MSG_ERROR([The memory header is required for peopt])])
AC_CHECK_HEADERS(functional,,
    [AC_MSG_ERROR([The functional header is required for peopt])])
AC_CHECK_HEADERS(algorithm,,
    [AC_MSG_ERROR([The algorithm header is required for peopt])])
AC_CHECK_HEADERS(numeric,,
    [AC_MSG_ERROR([The numeric header is required for peopt])])
AC_LANG_POP(C++)

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_ARG_WITH([headers-only],
    AS_HELP_STRING(
        [--with-headers-only@<:@=ARG@:>@],
        [Skip the build process and only install a base set of headers.
        When this flag is set, the flags below are ignored.  (ARG=yes).
        Do a full build (ARGS="no" || "")]
    ),
    [if test "x$withval"="xyes"; then headers_only="yes"; fi]
)

# Checks for BLAS.
m4_include([m4/ax_blas.m4])
if test "x${headers_only}" != "xyes"; then
    AX_BLAS(have_blas=true)
fi

# Checks for LAPACK
m4_include([m4/ax_lapack.m4])
if test "x${headers_only}" != "xyes" && test "x$have_blas" == "xtrue"; then
    AX_LAPACK(have_lapack=true)
fi
AM_CONDITIONAL(HAVE_LAPACK,test "x$have_lapack" = "xtrue")

# Figure out the Fortran demangling
if test "x$have_blas" == "xtrue"; then
    AC_F77_WRAPPERS
fi

# Checks for OpenMP
m4_include([m4/ax_openmp.m4])
if test "x$headers_only" != "xyes"; then
    AC_LANG_PUSH(C++)
    AX_OPENMP()
    AC_SUBST(OPENMP_CXXFLAGS)
    AC_LANG_POP(C++)
fi

# Checks for JSONCPP
AC_LANG_PUSH(C++)
AC_ARG_WITH(jsoncpp,
    AS_HELP_STRING(
        [--with-jsoncpp@<:@=ARG@:>@],
        [Use jsoncpp from a standard location (ARG=yes || "").
        Use jsoncpp from a specified location (ARG=<path>).
        Do not use jsoncpp (ARGS="no").]
    )
)
if test "x${headers_only}" != "xyes"; then
    # First, check for the header in the standard location
    if test "x$with_jsoncpp" == "xyes" || test "x$with_jsoncpp" == "x"; then
        AC_CHECK_HEADERS(json/json.h,,
            [AC_MSG_ERROR(
                [jsoncpp support enabled, but the header could not be found])])
        JSONCPP_LIBS="-ljson"
        jsoncpp_header=true

    # Then, check a custom location
    elif test "x$with_jsoncpp" != "xno"; then
        jsoncpp_path="$with_jsoncpp"
        AC_CHECK_HEADERS($jsoncpp_path/include/json/json.h,,
            [AC_MSG_ERROR(
                [jsoncpp support enabled, but the header could not be found])])
        JSONCPP_CPPFLAGS="-I$jsoncpp_path/include"
        AC_SUBST(JSONCPP_CPPFLAGS)
        JSONCPP_LDFLAGS="-L$jsoncpp_path/lib"
        AC_SUBST(JSONCPP_LDFLAGS)
        JSONCPP_LIBS="-ljson"
        AC_SUBST(JSONCPP_LIBS)
        jsoncpp_header=true
    fi
fi
# If we've succesfully found the header, check for the library
if test "x$jsoncpp_header" == "xtrue"
then
    old_CPPFLAGS=$CPPFLAGS 
    old_LDFLAGS=$LDFLAGS
    old_LIBS=$LIBS 
    CPPFLAGS="$CPPFLAGS $JSONCPP_CPPFLAGS"
    LDFLAGS="$LDFLAGS $JSONCPP_LDFLAGS"
    LIBS="$LIBS $JSONCPP_LIBS"
    AC_CACHE_CHECK(
        [whether the jsoncpp library is available],
        [jsoncpp_cv_lib],
        [AC_LINK_IFELSE(
            [AC_LANG_PROGRAM([#include <json/json.h>],[Json::Value dummy])],
            [jsoncpp_cv_lib=yes],                                           
            [AC_MSG_ERROR(
                [jsoncpp support enabled, but the library could not be found])])
        ]
    )
    CPPFLAGS=$old_CPPFLAGS 
    LDFLAGS=$old_LDFLAGS
    LIBS=$old_LIBS 
fi
AM_CONDITIONAL(HAVE_JSONCPP,test "x$jsoncpp_cv_lib" = "xyes")
AC_LANG_POP(C++)

# Checks for the Boost libraries
if test "x${headers_only}" != "xyes"; then
    m4_include([m4/ax_boost_base.m4])
    m4_include([m4/ax_boost_unit_test_framework.m4])
    AX_BOOST_BASE(,boost_base=true)
    if test "x$boost_base" = "xtrue"
    then
        AX_BOOST_UNIT_TEST_FRAMEWORK()
    fi
fi
AM_CONDITIONAL(HAVE_BOOST,test "x$BOOST_UNIT_TEST_FRAMEWORK_LIB" != "x")

AC_CONFIG_FILES([Makefile
                 src/peopt/Makefile
                 src/unit/Makefile
                 src/examples/rosenbrock/Makefile
                 src/examples/rosenbrock_easy/Makefile])
AC_OUTPUT
