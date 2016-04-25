# Add a module directory in case we need to find some packages
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
    "${CMAKE_SOURCE_DIR}/src/cmake/Modules/")
    
# Add some useful macros for examples and unit testing 
include(Optizelle)

# See if we're going to build the documentation for the package
set(ENABLE_DOCUMENTATION OFF CACHE BOOL "Enable documentation generation?")
if(ENABLE_DOCUMENTATION)
    add_subdirectory(doc)

# If we're not building the documentation, hide some of the parameters
else()
    set(ENABLE_A4_PAPER OFF CACHE BOOL "Generate A4 sized documentation?  If false, the document is letter sized." FORCE)
    mark_as_advanced(FORCE ENABLE_A4_PAPER)
endif()
        
# Figure out if we're doing a barebones header only installation
# or if we build the library.
set(ENABLE_CPP OFF CACHE BOOL "Enable libraries for C++?")
if(ENABLE_CPP)

    # Force C++14 on all targets
    set(CMAKE_CXX_STANDARD 14)
    set_property(GLOBAL PROPERTY CXX_STANDARD_REQUIRED)

    # Turn on position independent code.  This simplifies building both static
    # and shared libraries.
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)

    # Grab our third party library dependencies
    add_subdirectory(thirdparty)

    # Compile our C++ pieces 
    add_subdirectory(cpp)

    # Compile our Python pieces
    add_subdirectory(python)

    # Compile our MATLAB/Octave pieces
    add_subdirectory(matlab_octave)
    
    # Compile our unit tests
    add_subdirectory(unit)

    # Compile our examples 
    add_subdirectory(examples)

# Hide some options if the library is not enabled
else()
    set(ENABLE_BUILD_JSONCPP OFF CACHE BOOL "Build jsoncpp from source?" FORCE)
    set(ENABLE_BUILD_BLAS_AND_LAPACK OFF CACHE BOOL
        "Build BLAS and LAPACK from source?" FORCE)
    set(ENABLE_OPENMP OFF CACHE BOOL "Enable OpenMP?" FORCE)
    set(ENABLE_CPP_EXAMPLES OFF CACHE BOOL "Enable examples for C++?" FORCE)
    set(ENABLE_CPP_UNIT OFF CACHE BOOL "Enable unit tests for C++?" FORCE)
    set(ENABLE_PYTHON OFF CACHE BOOL "Enable the Python build?" FORCE)
    set(ENABLE_PYTHON_EXAMPLES OFF CACHE BOOL "Enable examples for Python?"
        FORCE)
    set(ENABLE_PYTHON_UNIT OFF CACHE BOOL "Enable unit tests for Python?" FORCE)
    set(ENABLE_MATLAB OFF CACHE BOOL "Enable the MATLAB build?" FORCE)
    set(ENABLE_MATLAB_EXAMPLES OFF CACHE BOOL "Enable examples for MATLAB?"
        FORCE)
    set(ENABLE_MATLAB_UNIT OFF CACHE BOOL "Enable unit tests for MATLAB?" FORCE)
    mark_as_advanced(FORCE
        ENABLE_BUILD_JSONCPP
        ENABLE_BUILD_BLAS_AND_LAPACK
        ENABLE_CPP_EXAMPLES
        ENABLE_CPP_UNIT
        ENABLE_PYTHON
        ENABLE_PYTHON_EXAMPLES
        ENABLE_PYTHON_UNIT
        ENABLE_MATLAB
        ENABLE_MATLAB_EXAMPLES
        ENABLE_MATLAB_UNIT
        ENABLE_OPENMP)
endif()
