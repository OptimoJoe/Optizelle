# This is the minimum required for the library install
#cmake_minimum_required(VERSION 2.8.9)

# This is the minimum required for the header only install
cmake_minimum_required(VERSION 2.4)
project(peopt)

# Add a module directory in case we need to find some packages
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
    "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# If we don't compile anything, install the basic headers
install(FILES
    ${CMAKE_SOURCE_DIR}/cpp/peopt/peopt.h
    ${CMAKE_SOURCE_DIR}/cpp/peopt/json.h
    ${CMAKE_SOURCE_DIR}/cpp/peopt/linalg.h
    DESTINATION include/peopt)

# See if we're going to build the documentation for the package
set(ENABLE_DOCUMENTATION OFF CACHE BOOL "Enable documentation generation?")
if(ENABLE_DOCUMENTATION)
    include(UseLATEX)
    add_subdirectory(doc)
endif()
        
# Figure out if we're doing a barebones header only installation
# or if we build the library.
set(ENABLE_CPP_LIBRARY OFF CACHE BOOL "Enable libraries for C++?")
if(ENABLE_CPP_LIBRARY)
    # Show the available options once the library is enabled
    mark_as_advanced(CLEAR
        ENABLE_CPP_EXAMPLES
        ENABLE_CPP_UNIT
        ENABLE_PYTHON
        ENABLE_MATLAB
        ENABLE_OPENMP)

    # Turn on position independent code immediately.  This simplifies
    # building both static and shared libraries.
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)

    # In order to build the library, especially the complete library with
    # all dependencies, find the required dependencies 
    enable_language(Fortran)
    find_package(JsonCpp REQUIRED)
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)

    # Figure out if we should look for OpenMP
    set(ENABLE_OPENMP OFF CACHE BOOL "Enable OpenMP?")
    if(ENABLE_OPENMP)
        # If we enable OpenMP, go ahead and make it required and find the
        # library.
        find_package(OpenMP REQUIRED)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
    endif()

    # Add in the directory with the library sources
    add_subdirectory(cpp/peopt)

    # Set the include path to include cpp since we will likely need the
    # peopt header file
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/cpp)

    # Check if we want to build the C++ examples
    set(ENABLE_CPP_EXAMPLES OFF CACHE BOOL "Enable examples for C++?")
    if(ENABLE_CPP_EXAMPLES)
        add_subdirectory(cpp/examples/diagnostic_checks)
        add_subdirectory(cpp/examples/rosenbrock_custom_VS)
        add_subdirectory(cpp/examples/rosenbrock)
        add_subdirectory(cpp/examples/quadratic)
        add_subdirectory(cpp/examples/simple_inequality)
        add_subdirectory(cpp/examples/simple_infeasible_inequality)
        add_subdirectory(cpp/examples/simple_quadratic_cone)
        add_subdirectory(cpp/examples/simple_sdp_cone)
        add_subdirectory(cpp/examples/simple_equality)
        add_subdirectory(cpp/examples/nw_sqp_exercise)
        add_subdirectory(cpp/examples/simple_constrained)
        add_subdirectory(cpp/examples/sdpa_sparse_format)
        add_subdirectory(cpp/examples/preconditioned_nonlinear_cg)
    endif()
    
    # Check if we want to build the unit tests
    set(ENABLE_CPP_UNIT OFF CACHE BOOL "Enable unit tests for C++?")
    if(ENABLE_CPP_UNIT)
        # Assuming we enable the unit tests, we require the
        # boost_unit_test_framework library and headers
        find_package(BoostUnitTest REQUIRED)
        add_subdirectory(cpp/unit)
    endif()
    
    # Check if we want to build the Python libraries
    set(ENABLE_PYTHON OFF CACHE BOOL "Enable the Python build?")
    if(ENABLE_PYTHON)
        # Make sure we can find the Python libraries and then add the
        # directory with the sources.
        find_package(PythonLibs)
        add_subdirectory(python/peopt)
        
        # Check if we want to add examples for Python
        set(ENABLE_PYTHON_EXAMPLES OFF CACHE BOOL "Enable examples for Python?")
        mark_as_advanced(CLEAR ENABLE_PYTHON_EXAMPLES)
        if(ENABLE_PYTHON_EXAMPLES)
            add_subdirectory(python/examples/diagnostic_checks)
            add_subdirectory(python/examples/rosenbrock)
            add_subdirectory(python/examples/simple_inequality)
        endif()
    else()
        set(ENABLE_PYTHON_EXAMPLES OFF CACHE BOOL "Enable examples for Python?"
            FORCE)
        mark_as_advanced(FORCE ENABLE_PYTHON_EXAMPLES)
    endif()
    
    # Check if we want to build the MATLAB interface 
    set(ENABLE_MATLAB OFF CACHE BOOL "Enable the MATLAB build?")
    if(ENABLE_MATLAB)
        # Make sure we can find the MATLAB libraries and then add the
        # directory with the sources.
        find_package(MATLAB)
        add_subdirectory(matlab/peopt)

        # Check if we want to add examples for MATLAB 
        set(ENABLE_MATLAB_EXAMPLES OFF CACHE BOOL "Enable examples for MATLAB?")
        mark_as_advanced(CLEAR ENABLE_MATLAB_EXAMPLES)
        if(ENABLE_MATLAB_EXAMPLES)
            add_subdirectory(matlab/examples/diagnostic_checks)
            add_subdirectory(matlab/examples/rosenbrock)
            add_subdirectory(matlab/examples/rosenbrock_ineq)
            add_subdirectory(matlab/examples/parest)
        endif()
    else()
        set(ENABLE_MATLAB_EXAMPLES OFF CACHE BOOL "Enable examples for MATLAB?"
            FORCE)
        mark_as_advanced(FORCE ENABLE_MATLAB_EXAMPLES)
    endif()

    # Build the library
    if(ENABLE_PYTHON)
        add_library(peopt_shared SHARED
            $<TARGET_OBJECTS:peopt_cpp>
            $<TARGET_OBJECTS:peopt_python>)
        add_library(peopt_shared_complete SHARED
            $<TARGET_OBJECTS:peopt_cpp>
            $<TARGET_OBJECTS:peopt_python>)
        add_library(peopt_static STATIC 
            $<TARGET_OBJECTS:peopt_cpp>
            $<TARGET_OBJECTS:peopt_python>)
    else()
        add_library(peopt_shared SHARED
            $<TARGET_OBJECTS:peopt_cpp>)
        add_library(peopt_shared_complete SHARED
            $<TARGET_OBJECTS:peopt_cpp>)
        add_library(peopt_static STATIC 
            $<TARGET_OBJECTS:peopt_cpp>)
    endif()
    target_link_libraries(peopt_shared_complete
        ${JSONCPP_LIBRARIES}
        ${LAPACK_LIBRARIES}
        ${BLAS_LIBRARIES})

    # This is required for Windows to see an import library, which is
    # required to build a dll.
    target_link_libraries(peopt_shared peopt_shared_complete)

    set_target_properties(peopt_shared peopt_static
        PROPERTIES OUTPUT_NAME peopt)
    set_target_properties(peopt_shared_complete 
        PROPERTIES OUTPUT_NAME peopt_complete)
    install(TARGETS peopt_shared peopt_shared_complete peopt_static
        DESTINATION lib)

# Hide some options if the library is not enabled
else()
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
        ENABLE_CPP_EXAMPLES
        ENABLE_CPP_UNIT
        ENABLE_PYTHON
        ENABLE_PYTHON_EXAMPLES
        ENABLE_MATLAB
        ENABLE_MATLAB_EXAMPLES
        ENABLE_OPENMP)
endif()
