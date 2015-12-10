# Set the project name
project(Optizelle)

# This is the minimum required for the install.  Basically, we need
# CXX_STANDARD_REQUIRED and this looks to have been added at 3.1.  I think
# 3.1.3 is what the current bug fix is, so that looks good enough. 
cmake_minimum_required(VERSION 3.1.3 FATAL_ERROR)

# Set the CPack options
set(CPACK_GENERATOR "TGZ")
set(CPACK_PACKAGE_VENDOR "OptimoJoe")
set(CPACK_PACKAGE_NAME "Optizelle")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "1")
set(CPACK_PACKAGE_VERSION_PATCH "5")
set(CPACK_SOURCE_IGNORE_FILES "${CMAKE_BINARY_DIR}" "${CMAKE_SOURCE_DIR}/.git")
include(CPack)
    
# Enable ctest
enable_testing()

# Add the source directories for all of the ports
add_subdirectory(src)
