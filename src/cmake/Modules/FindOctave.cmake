# Handle pkg-config
find_package(PkgConfig)
pkg_check_modules(PC_OCTAVE QUIET octave)
set(OCTAVE_DEFINITIONS ${PC_OCTAVE_CFLAGS_OTHER})

# Find the header file
find_path(OCTAVE_INCLUDE_DIR mex.h
    HINTS ${PC_OCTAVE_INCLUDEDIR} ${PC_OCTAVE_INCLUDE_DIRS}
    PATH_SUFFIXES matlab )
set(OCTAVE_INCLUDE_DIRS ${OCTAVE_INCLUDE_DIR} )

# Find the library
find_library(OCTAVE_LIBRARY
    NAMES octinterp liboctinterp
    HINTS ${PC_OCTAVE_LIBDIR} ${PC_OCTAVE_LIBRARY_DIRS} )
set(OCTAVE_LIBRARIES ${OCTAVE_LIBRARY} )

# Handle the quietly and required arguments for the library and include
# directory.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OCTAVE DEFAULT_MSG
    OCTAVE_LIBRARY OCTAVE_INCLUDE_DIR)
mark_as_advanced(OCTAVE_INCLUDE_DIR OCTAVE_LIBRARY)

# Find the Octave executable
set(OCTAVE_EXECUTABLE_UNDEFINED "OCTAVE_EXECUTABLE-UNDEFINED")
set(OCTAVE_EXECUTABLE ${OCTAVE_EXECUTABLE_UNDEFINED} CACHE FILEPATH
    "Octave executable with full path.")
# Make sure that it is defined and throw an error message if it is not
if(OCTAVE_EXECUTABLE STREQUAL OCTAVE_EXECUTABLE_UNDEFINED)
    mark_as_advanced(CLEAR OCTAVE_EXECUTABLE)
    message(FATAL_ERROR
        "The Octave executable with full path must be defined in the variable OCTAVE_EXECUTABLE.")
endif()
mark_as_advanced(FORCE OCTAVE_EXECUTABLE)
