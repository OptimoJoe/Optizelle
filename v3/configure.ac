# Initialization
AC_PREREQ(2.59)
AC_INIT(peopt, 0.3.0, josyoun@sandia.gov)
AC_CONFIG_AUX_DIR(build-aux)
AM_INIT_AUTOMAKE(peopt, 0.3.0, josyoun@sandia.gov)

# Checks for programs.
AC_PROG_CXX       # Needed for the primary program
AC_PROG_CC        # Needed for the tests below.  For example, OpenMP
AC_PROG_F77       # Needed for testing BLAS and LAPACK
AC_PROG_FC        # Needed for name mangling
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

# Figure out the Fortran demangling
AC_FC_WRAPPERS

# Checks for BLAS and LAPACK
m4_include([m4/ax_blas.m4])
m4_include([m4/ax_lapack.m4])
AX_BLAS()
AX_LAPACK(have_lapack=true)
AM_CONDITIONAL(HAVE_LAPACK,test "$have_lapack" = "true")

# Checks for OpenMP
m4_include([m4/ax_openmp.m4])
AC_LANG_PUSH(C++)
AX_OPENMP()
AC_SUBST(OPENMP_CXXFLAGS)
AC_LANG_POP(C++)

# Checks for JSONCPP
AC_LANG_PUSH(C++)
AC_ARG_WITH(jsoncpp,[
  --with-jsoncpp@<:@=ARG@:>@    use jsoncpp from a standard location (ARG=yes)
                          use jsoncpp from a specified location (ARG=<path>)
                          do not use jsoncpp (ARG="")])
# First, check for the header
if test "x$with_jsoncpp" == "xyes"
then
    AC_CHECK_HEADERS(json/json.h,,
        [AC_MSG_ERROR(
            [jsoncpp support enabled, but the header could not be found])])
    JSONCPP_LIBS="-ljson"
elif test "x$with_jsoncpp" != "x"
then
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
fi
# If we've succesfully found the header, check for the library
if test "x$with_jsoncpp" != "x"
then
    old_CPPFLAGS=$CPPFLAGS 
    old_LDFLAGS=$LDFLAGS
    old_LIBS=$LIBS 
    CPPFLAGS="$CPPFLAGS $JSONCPP_CPPFLAGS"
    LDFLAGS="$LDFLAGS $JSONCPP_LDFLAGS"
    LIBS="$LIBS $JSONCPP_LIBS"
    AC_CACHE_CHECK(
        [whether the jsoncpp library is available],
        [jsoncpp_lib],
        [AC_LINK_IFELSE(
            [AC_LANG_PROGRAM([#include <json/json.h>],[Json::Value dummy])],
            [jsoncpp_lib=yes],                                           
            [AC_MSG_ERROR(
                [jsoncpp support enabled, but the library could not be found])])
        ]
    )
    CPPFLAGS=$old_CPPFLAGS 
    LDFLAGS=$old_LDFLAGS
    LIBS=$old_LIBS 
fi
AM_CONDITIONAL(HAVE_JSONCPP,test "x$jsoncpp_lib" = "xyes")
AC_LANG_POP(C++)

# Checks for the Boost libraries
m4_include([m4/ax_boost_base.m4])
m4_include([m4/ax_boost_unit_test_framework.m4])
AX_BOOST_BASE(,boost_base=true)
if test "x$boost_base" = "xtrue"
then
    AX_BOOST_UNIT_TEST_FRAMEWORK()
fi
AM_CONDITIONAL(HAVE_BOOST,test "x$BOOST_UNIT_TEST_FRAMEWORK_LIB" != "x")

AC_CONFIG_FILES([Makefile
                 src/peopt/Makefile
                 src/unit/Makefile
                 src/examples/rosenbrock/Makefile
                 src/examples/rosenbrock_easy/Makefile])
AC_OUTPUT
