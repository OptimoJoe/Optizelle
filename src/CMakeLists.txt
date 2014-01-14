# Add a module directory in case we need to find some packages
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
    "${CMAKE_SOURCE_DIR}/src/cmake/Modules/")

# See if we're going to build the documentation for the package
set(ENABLE_DOCUMENTATION OFF CACHE BOOL "Enable documentation generation?")
if(ENABLE_DOCUMENTATION)
    add_subdirectory(doc)
endif()
        
# Figure out if we're doing a barebones header only installation
# or if we build the library.
set(ENABLE_CPP_LIBRARY OFF CACHE BOOL "Enable libraries for C++?")
if(ENABLE_CPP_LIBRARY)

    # Add some useful macros for unit testing
    include(AddOptizelleTestCpp)
    include(AddOptizelleSDPACpp)
    include(AddOptizelleExampleCpp)
    include(AddOptizelleExampleSupporting)

    # Figure out if we should look for OpenMP
    mark_as_advanced(CLEAR ENABLE_OPENMP)
    set(ENABLE_OPENMP OFF CACHE BOOL "Enable OpenMP?")
    if(ENABLE_OPENMP)
        # If we enable OpenMP, go ahead and make it required and find the
        # library.
        find_package(OpenMP REQUIRED)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
    endif()

    # Set the Optizelle include directories
    set(OPTIZELLE_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/cpp)

    # Grab our third party library dependencies
    add_subdirectory(thirdparty)

    # Compile our C++ pieces 
    add_subdirectory(cpp)

    # Compile our Python pieces
    add_subdirectory(python)

    # Compile our MATLAB pieces
    add_subdirectory(matlab)

# Hide some options if the library is not enabled
else()
    set(ENABLE_BUILD_JSONCPP OFF CACHE BOOL "Build jsoncpp from source?" FORCE)
    set(ENABLE_BUILD_BLAS_AND_LAPACK OFF CACHE BOOL
        "Build BLAS and LAPACK from source?" FORCE)
    set(ENABLE_CPP_EXAMPLES OFF CACHE BOOL "Enable examples for C++?" FORCE)
    set(ENABLE_CPP_UNIT OFF CACHE BOOL "Enable unit tests for C++?" FORCE)
    set(ENABLE_PYTHON OFF CACHE BOOL "Enable the Python build?" FORCE)
    set(ENABLE_PYTHON_EXAMPLES OFF CACHE BOOL "Enable examples for Python?"
        FORCE)
    set(ENABLE_MATLAB OFF CACHE BOOL "Enable the MATLAB build?" FORCE)
    set(ENABLE_MATLAB_EXAMPLES OFF CACHE BOOL "Enable examples for MATLAB?"
        FORCE)
    set(ENABLE_OPENMP OFF CACHE BOOL "Enable OpenMP?" FORCE)
    mark_as_advanced(FORCE
        ENABLE_BUILD_JSONCPP
        ENABLE_BUILD_BLAS_AND_LAPACK
        ENABLE_CPP_EXAMPLES
        ENABLE_CPP_UNIT
        ENABLE_PYTHON
        ENABLE_PYTHON_EXAMPLES
        ENABLE_MATLAB
        ENABLE_MATLAB_EXAMPLES
        ENABLE_OPENMP)
endif()
