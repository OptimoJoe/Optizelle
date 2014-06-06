# Set the project name
project(Optizelle)

# This is the minimum required for the install.  Basically, we need
# set(CMAKE_POSITION_INDEPENDENT_CODE ON) and this was added in 2.8.9.
cmake_minimum_required(VERSION 2.8.9)

# Set the CPack options
set(CPACK_GENERATOR "TGZ")
set(CPACK_PACKAGE_VENDOR "OptimoJoe")
set(CPACK_PACKAGE_NAME "Optizelle")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_SOURCE_IGNORE_FILES "${CMAKE_BINARY_DIR}" "${CMAKE_SOURCE_DIR}/.git")
include(CPack)
    
# Enable ctest
enable_testing()

# Add the source directories for all of the ports
add_subdirectory(src)
