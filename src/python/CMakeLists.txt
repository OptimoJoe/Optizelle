# Check if we want to build the Python libraries
mark_as_advanced(CLEAR ENABLE_PYTHON)
set(ENABLE_PYTHON OFF CACHE BOOL "Enable the Python build?")
if(ENABLE_PYTHON)
    # Common includes
    include_directories(${JSONCPP_INCLUDE_DIRS})
    include_directories(${PEOPT_INCLUDE_DIRS})

    # Make sure we can find the Python libraries and then add the
    # directory with the sources.
    find_package(PythonLibs)
    add_subdirectory(peopt)
    
    # Check if we want to add examples for Python
    set(ENABLE_PYTHON_EXAMPLES OFF CACHE BOOL "Enable examples for Python?")
    mark_as_advanced(CLEAR ENABLE_PYTHON_EXAMPLES)
    if(ENABLE_PYTHON_EXAMPLES)
        add_subdirectory(examples)
    endif()
else()
    set(ENABLE_PYTHON_EXAMPLES OFF CACHE BOOL "Enable examples for Python?"
        FORCE)
    mark_as_advanced(FORCE ENABLE_PYTHON_EXAMPLES)
endif()
