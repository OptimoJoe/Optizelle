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
set(CPACK_PACKAGE_LONG_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(CPACK_SOURCE_IGNORE_FILES "${CMAKE_BINARY_DIR}" "${CMAKE_SOURCE_DIR}/.git")
set(CPACK_PACKAGE_CONTACT "Joseph Young <joe@optimojoe.com>")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Optizelle [op-tuh-zel] is an open source software library designed to solve general purpose nonlinear optimization problems")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_BINARY_DIR}/licenses/LICENSE.txt")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "Optizelle")
if(WIN32)
    set(CPACK_GENERATOR "WIX")
    set(CPACK_WIX_UPGRADE_GUID "9FB8001A-10E0-4FD1-9DC7-4B5415D6676E")
    set(CPACK_WIX_PATCH_FILE "${CMAKE_SOURCE_DIR}/src/cmake/WiX/patch.xml")
    set(CPACK_WIX_EXTENSIONS "WixUtilExtension")
elseif(APPLE)
    # Create a dmg file
    set(CPACK_GENERATOR "DragNDrop")

    # Set the installation folder
    set(CPACK_PACKAGING_INSTALL_PREFIX "/${CPACK_PACKAGE_LONG_NAME}")

    # Make sure to enable this when building a package
    set(CPACK_DRAGNDROP OFF CACHE BOOL "Building a DragNDrop package?")
    if(CPACK_DRAGNDROP)
        # Setup the readme file
        set(README_NAME "src/cmake/DragNDrop/ReadMe.txt")
        configure_file(
            "${CMAKE_SOURCE_DIR}/${README_NAME}"
            "${CMAKE_BINARY_DIR}/${README_NAME}")
        install(FILES "${CMAKE_BINARY_DIR}/${README_NAME}" DESTINATION ..)

        # Setup the plist file
        set(PLIST_NAME "src/cmake/DragNDrop/com.optimojoe.optizelle.plist")
        configure_file(
            "${CMAKE_CURRENT_SOURCE_DIR}/${PLIST_NAME}"
            "${CMAKE_CURRENT_BINARY_DIR}/${PLIST_NAME}")

        # Correct an absurd number of hard coded paths
        macro(fix_macos_lib name lib)
            add_custom_target(fix_${name}_internal_${lib}
                ALL
                DEPENDS LAPACK
                COMMAND install_name_tool -change /opt/local/lib/libgcc/lib${lib}.dylib /Applications/${CPACK_PACKAGE_LONG_NAME}/lib/lib${lib}.dylib ${CMAKE_BINARY_DIR}/thirdparty/lib/lib${name}.dylib
                COMMENT "Fixing the path for lib${lib}.dylib inside of lib${name}.dylib")
        endmacro()
        macro(fix_macos_lib_alt name lib)
            add_custom_target(fix_${name}_internal_${lib}
                ALL
                DEPENDS LAPACK
                COMMAND install_name_tool -change /opt/local/lib/gcc43/lib${lib}.dylib /Applications/${CPACK_PACKAGE_LONG_NAME}/lib/lib${lib}.dylib ${CMAKE_BINARY_DIR}/thirdparty/lib/lib${name}.dylib
                COMMENT "Fixing the path for lib${lib}.dylib inside of lib${name}.dylib")
        endmacro()
        fix_macos_lib(gfortran.3 quadmath.0)
        fix_macos_lib(gfortran.3 gcc_s.1)
        fix_macos_lib(quadmath.0 gcc_s.1)
        fix_macos_lib(blas.3.6.1 gfortran.3)
        fix_macos_lib_alt(blas.3.6.1 gcc_s.1)
        fix_macos_lib(lapack.3.6.1 gfortran.3)
        fix_macos_lib_alt(lapack.3.6.1 gcc_s.1)
    endif()
elseif(UNIX)
    set(CPACK_GENERATOR "TGZ" "RPM" "DEB")
endif()
include(CPack)
    
# Enable ctest
enable_testing()

# Add the source directories for all of the ports
add_subdirectory(src)
