# Handle pkg-config
find_package(PkgConfig)
pkg_check_modules(PC_JSONCPP QUIET jsoncpp)
set(JSONCPP_DEFINITIONS ${PC_JSONCPP_CFLAGS_OTHER})

# Find the header file
find_path(JSONCPP_INCLUDE_DIR json/json.h
    HINTS ${PC_JSONCPP_INCLUDEDIR} ${PC_JSONCPP_INCLUDE_DIRS}
    PATH_SUFFIXES json )
set(JSONCPP_INCLUDE_DIRS ${JSONCPP_INCLUDE_DIR} )

# Find the library
find_library(JSONCPP_LIBRARY
    NAMES jsoncpp json libjson
    HINTS ${PC_JSONCPP_LIBDIR} ${PC_JSONCPP_LIBRARY_DIRS} )
set(JSONCPP_LIBRARIES ${JSONCPP_LIBRARY} )

# Handle the quietly and required arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(JsonCpp DEFAULT_MSG
    JSONCPP_LIBRARY JSONCPP_INCLUDE_DIR)
mark_as_advanced(JSONCPP_INCLUDE_DIR JSONCPP_LIBRARY)
