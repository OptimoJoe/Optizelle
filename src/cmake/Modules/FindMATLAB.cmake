# Handle pkg-config 
find_package(PkgConfig)
pkg_check_modules(PC_MATLAB QUIET matlab)
set(MATLAB_DEFINITIONS ${PC_MATLAB_CFLAGS_OTHER})

# Find the header file
find_path(MATLAB_INCLUDE_DIR mex.h
    HINTS ${PC_MATLAB_INCLUDEDIR} ${PC_MATLAB_INCLUDE_DIRS}
    PATH_SUFFIXES matlab )
set(MATLAB_INCLUDE_DIRS ${MATLAB_INCLUDE_DIR} )

# Find the library
find_library(MATLAB_LIBRARY
    NAMES mex libmex
    HINTS ${PC_MATLAB_LIBDIR} ${PC_MATLAB_LIBRARY_DIRS} )
set(MATLAB_LIBRARIES ${MATLAB_LIBRARY} )

# Handle the quietly and required arguments for the library and include
# directory.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MATLAB DEFAULT_MSG
    MATLAB_LIBRARY MATLAB_INCLUDE_DIR) 
mark_as_advanced(MATLAB_INCLUDE_DIR MATLAB_LIBRARY)

# Find the mex file extension
set(MEX_UNDEFINED "MATLAB_MEX_EXTENSION-UNDEFINED")
set(MATLAB_MEX_EXTENSION ${MEX_UNDEFINED} CACHE STRING
    "MATLAB mex file extension.  Type 'mexext' in MATLAB.")
# Make sure that it is defined and throw an error message if it is not
if(MATLAB_MEX_EXTENSION STREQUAL MEX_UNDEFINED) 
    mark_as_advanced(CLEAR MATLAB_MEX_EXTENSION)
    message(FATAL_ERROR
        "The MATLAB mex extension must be defined in the variable MATLAB_MEX_EXTENSION.  Type 'mexext' in MATLAB.")
endif()
mark_as_advanced(FORCE MATLAB_MEX_EXTENSION)

# Find the MATLAB executable 
set(MEX_UNDEFINED "MATLAB_EXECUTABLE-UNDEFINED")
set(MATLAB_EXECUTABLE ${MEX_UNDEFINED} CACHE STRING
    "MATLAB executable with full path.")
# Make sure that it is defined and throw an error message if it is not
if(MATLAB_EXECUTABLE STREQUAL MEX_UNDEFINED) 
    mark_as_advanced(CLEAR MATLAB_EXECUTABLE)
    message(FATAL_ERROR
        "The MATLAB executable with full path must be defined in the variable MATLAB_EXECUTABLE.")
endif()
mark_as_advanced(FORCE MATLAB_EXECUTABLE)

# Find the MATLAB executable 
set(MEX_UNDEFINED "MATLAB_RUN_FLAG-UNDEFINED")
set(MATLAB_RUN_FLAG ${MEX_UNDEFINED} CACHE STRING
    "MATLAB command line flag to run a script (-r in MATLAB and --eval in Octave).")
# Make sure that it is defined and throw an error message if it is not
if(MATLAB_RUN_FLAG STREQUAL MEX_UNDEFINED) 
    mark_as_advanced(CLEAR MATLAB_RUN_FLAG)
    message(FATAL_ERROR
        "The MATLAB command line flag to run a script (-r in MATLAB and --eval in Octave) must be defined in the variable MATLAB_RUN_FLAG.")
endif()
mark_as_advanced(FORCE MATLAB_RUN_FLAG)
