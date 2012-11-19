# Handle pkg-config 
find_package(PkgConfig)
pkg_check_modules(PC_BLAS QUIET blas)
set(BLAS_DEFINITIONS ${PC_BLAS_CFLAGS_OTHER})

# Find the library
find_library(BLAS_LIBRARY
    NAMES blas libblas
    HINTS ${PC_BLAS_LIBDIR} ${PC_BLAS_LIBRARY_DIRS} )
set(BLAS_LIBRARIES ${BLAS_LIBRARY} )

# Handle the quietly and required arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BLAS DEFAULT_MSG BLAS_LIBRARY)
mark_as_advanced(BLAS_LIBRARY)
