# Enable C++
enable_language(CXX)

# Set common includes
include_directories(${OPTIZELLE_INCLUDE_DIRS})
include_directories(${JSONCPP_INCLUDE_DIRS})

# Build the Optizelle library 
add_subdirectory(optizelle)

# Check if we want to build the unit tests
mark_as_advanced(CLEAR ENABLE_CPP_UNIT)
set(ENABLE_CPP_UNIT OFF CACHE BOOL "Enable unit tests for C++?")
if(ENABLE_CPP_UNIT)
    # Basic unit tests and utility functions
    add_subdirectory(unit)

    # Unit tests based on examples
    add_subdirectory(examples)
endif()

# Check if we want to build the C++ examples.  
mark_as_advanced(CLEAR ENABLE_CPP_EXAMPLES)
set(ENABLE_CPP_EXAMPLES OFF CACHE BOOL "Enable examples for C++?")
if(ENABLE_CPP_EXAMPLES AND NOT ENABLE_CPP_UNIT)
    add_subdirectory(examples)
endif()
