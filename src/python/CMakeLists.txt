# Check if we want to build the Python libraries
mark_as_advanced(CLEAR ENABLE_PYTHON)
set(ENABLE_PYTHON OFF CACHE BOOL "Enable the Python build?")
if(ENABLE_PYTHON)
    # Common includes
    include_directories(${JSONCPP_INCLUDE_DIRS})
    include_directories(${OPTIZELLE_INCLUDE_DIRS})

    # Make sure we can find the Python libraries and then add the
    # directory with the sources.
    #set(Python_ADDITIONAL_VERSIONS 2.7)
    find_package(PythonLibs REQUIRED)
    find_package(PythonInterp REQUIRED)
    add_subdirectory(Optizelle)
endif()
