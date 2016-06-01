# Set the project name
project(Optizelle)

# This is the minimum required for the install.  Basically, we need
# CXX_STANDARD_REQUIRED and this looks to have been added at 3.1.  I think
# 3.1.3 is what the current bug fix is, so that looks good enough. 
cmake_minimum_required(VERSION 3.1.3 FATAL_ERROR)

# Put together the complete license
add_subdirectory(licenses)

# Install the authors file 
install(FILES AUTHORS.txt DESTINATION share/optizelle)

# Set the CPack options
set(CPACK_GENERATOR "TGZ")
set(CPACK_PACKAGE_VENDOR "OptimoJoe")
set(CPACK_PACKAGE_NAME "Optizelle")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "1")
set(CPACK_PACKAGE_VERSION_PATCH "5")
set(CPACK_SOURCE_IGNORE_FILES "${CMAKE_BINARY_DIR}" "${CMAKE_SOURCE_DIR}/.git")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "Optizelle")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_BINARY_DIR}/licenses/LICENSE.txt")
set(CPACK_WIX_UPGRADE_GUID "9FB8001A-10E0-4FD1-9DC7-4B5415D6676E")
set(CPACK_WIX_PATCH_FILE "${CMAKE_SOURCE_DIR}/src/CMake/WiX/patch.xml")
set(CPACK_WIX_EXTENSIONS "WixUtilExtension")
include(CPack)
    
# Enable ctest
enable_testing()

# Add the source directories for all of the ports
add_subdirectory(src)
