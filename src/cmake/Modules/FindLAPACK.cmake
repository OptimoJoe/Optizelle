# Handle pkg-config
find_package(PkgConfig)
pkg_check_modules(PC_LAPACK QUIET lapack)
set(LAPACK_DEFINITIONS ${PC_LAPACK_CFLAGS_OTHER})

# Find the library
find_library(LAPACK_LIBRARY
    NAMES lapack liblapack
    HINTS ${PC_LAPACK_LIBDIR} ${PC_LAPACK_LIBRARY_DIRS} )
set(LAPACK_LIBRARIES ${LAPACK_LIBRARY} )

# Handle the quietly and required arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LAPACK DEFAULT_MSG LAPACK_LIBRARY)
mark_as_advanced(LAPACK_LIBRARY)
