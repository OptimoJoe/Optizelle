# Handle pkg-config 
find_package(PkgConfig)
pkg_check_modules(PC_BOOST_UNIT_TEST QUIET boost)
set(BOOST_UNIT_TEST_DEFINITIONS ${PC_BOOST_UNIT_TEST_CFLAGS_OTHER})

# Find the header file
find_path(BOOST_UNIT_TEST_INCLUDE_DIR boost/test/unit_test.hpp
    HINTS ${PC_BOOST_UNIT_TEST_INCLUDEDIR} ${PC_BOOST_UNIT_TEST_INCLUDE_DIRS}
    PATH_SUFFIXES boost/test)
set(BOOST_UNIT_TEST_INCLUDE_DIRS ${BOOST_UNIT_TEST_INCLUDE_DIR})

# Find the library
find_library(BOOST_UNIT_TEST_LIBRARY
    NAMES boost_unit_test_framework libboost_unit_test_framework
    HINTS ${PC_BOOST_UNIT_TEST_LIBDIR} ${PC_BOOST_UNIT_TEST_LIBRARY_DIRS})
set(BOOST_UNIT_TEST_LIBRARIES ${BOOST_UNIT_TEST_LIBRARY})

# Handle the quietly and required arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BOOST_UNIT_TEST DEFAULT_MSG
    BOOST_UNIT_TEST_LIBRARY BOOST_UNIT_TEST_INCLUDE_DIR)
mark_as_advanced(BOOST_UNIT_TEST_INCLUDE_DIR BOOST_UNIT_TEST_LIBRARY)
