# Common includes
include_directories(utility)

# Check if we want to build the unit tests
mark_as_advanced(CLEAR ENABLE_CPP_UNIT)
set(ENABLE_CPP_UNIT OFF CACHE BOOL "Enable unit tests for C++?")

# Basic unit tests and utility functions
add_subdirectory(capture_and_release)
add_subdirectory(derivative_check)
add_subdirectory(linear_algebra)
add_subdirectory(utility)

