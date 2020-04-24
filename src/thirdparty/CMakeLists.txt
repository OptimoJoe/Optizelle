project(third_party_libraries)

# Build jsoncpp
mark_as_advanced(CLEAR ENABLE_BUILD_JSONCPP)
set(ENABLE_BUILD_JSONCPP OFF CACHE BOOL "Build jsoncpp from source?")
if(ENABLE_BUILD_JSONCPP)
    # Determine the location of the JsonCpp archive
    set(JSONCPP_ARCHIVE_UNDEFINED "JSONCPP_ARCHIVE-UNDEFINED")
    set(JSONCPP_ARCHIVE ${JSONCPP_ARCHIVE_UNDEFINED} CACHE FILEPATH
        "JsonCpp archive with full path.")
    # Make sure that it is defined and throw an error message if it is not
    if(JSONCPP_ARCHIVE STREQUAL JSONCPP_ARCHIVE_UNDEFINED)
        mark_as_advanced(CLEAR JSONCPP_ARCHIVE)
        message(FATAL_ERROR
            "The JsonCpp archive with full path must be defined in the variable JSONCPP_ARCHIVE.")
    endif()
    mark_as_advanced(FORCE JSONCPP_ARCHIVE)

    include(ExternalProject)
    ExternalProject_Add(
        JsonCpp
        URL ${JSONCPP_ARCHIVE}
        CMAKE_ARGS
            -DCMAKE_BUILD_TYPE:STRING=RELEASE
            -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_BINARY_DIR}/thirdparty
            -DCMAKE_INSTALL_LIBDIR:STRING=lib
            -DCMAKE_AR:FILEPATH=${CMAKE_AR}
            -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
            -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
            -DCMAKE_LINKER:FILEPATH=${CMAKE_LINKER}
            -DCMAKE_MAKE_PROGRAM:FILEPATH=${CMAKE_MAKE_PROGRAM}
            -DCMAKE_NM:FILEPATH=${CMAKE_NM}
            -DCMAKE_OBJCOPY:FILEPATH=${CMAKE_OBJCOPY}
            -DCMAKE_OBJDUMP:FILEPATH=${CMAKE_OBJDUMP}
            -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
            -DCMAKE_STRIP:FILEPATH=${CMAKE_STRIP}
            -DBUILD_SHARED_LIBS:BOOL=ON
    )

    # Set the location of jsoncpp
    set(JSONCPP_LIBRARY
        "${CMAKE_BINARY_DIR}/thirdparty/lib/${CMAKE_SHARED_LIBRARY_PREFIX}jsoncpp${CMAKE_SHARED_LIBRARY_SUFFIX}"
        CACHE FILEPATH "Path to a library." FORCE)
    set(JSONCPP_LIBRARIES ${JSONCPP_LIBRARY})
    set(JSONCPP_INCLUDE_DIR
        "${CMAKE_BINARY_DIR}/thirdparty/include"
        CACHE
        PATH "Path to a file." FORCE)
    set(JSONCPP_INCLUDE_DIRS ${JSONCPP_INCLUDE_DIR})
    mark_as_advanced(JSONCPP_INCLUDE_DIR JSONCPP_LIBRARY)
else()
    find_package(JsonCpp REQUIRED)
endif()

# Pass all the library variables to the parent
set(JSONCPP_LIBRARY ${JSONCPP_LIBRARY} PARENT_SCOPE)
set(JSONCPP_LIBRARIES ${JSONCPP_LIBRARIES} PARENT_SCOPE)
set(JSONCPP_INCLUDE_DIR ${JSONCPP_INCLUDE_DIR} PARENT_SCOPE)
set(JSONCPP_INCLUDE_DIRS ${JSONCPP_INCLUDE_DIRS} PARENT_SCOPE)

# Build BLAS and LAPACK
mark_as_advanced(CLEAR ENABLE_BUILD_BLAS_AND_LAPACK)
set(ENABLE_BUILD_BLAS_AND_LAPACK OFF CACHE BOOL
    "Build BLAS and LAPACK from source?")
if(ENABLE_BUILD_BLAS_AND_LAPACK)
    # Determine the location of the LAPACK archive
    set(LAPACK_ARCHIVE_UNDEFINED "LAPACK_ARCHIVE-UNDEFINED")
    set(LAPACK_ARCHIVE ${LAPACK_ARCHIVE_UNDEFINED} CACHE FILEPATH
        "LAPACK archive with full path.")
    # Make sure that it is defined and throw an error message if it is not
    if(LAPACK_ARCHIVE STREQUAL LAPACK_ARCHIVE_UNDEFINED)
        mark_as_advanced(CLEAR LAPACK_ARCHIVE)
        message(FATAL_ERROR
            "The LAPACK archive with full path must be defined in the variable LAPACK_ARCHIVE.")
    endif()
    mark_as_advanced(FORCE LAPACK_ARCHIVE)

    include(ExternalProject)
    enable_language(Fortran)
    ExternalProject_Add(
        LAPACK
        URL ${LAPACK_ARCHIVE}
        CMAKE_ARGS
            -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_BINARY_DIR}/thirdparty
            -DBUILD_TESTING:BOOL=OFF
            -DCMAKE_INSTALL_LIBDIR:STRING=lib
            -DCMAKE_BUILD_TYPE:STRING=RELEASE
            -DBUILD_SHARED_LIBS:BOOL=ON
            -DCMAKE_AR:FILEPATH=${CMAKE_AR}
            -DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}
            -DCMAKE_Fortran_FLAGS:STRING=${CMAKE_Fortran_FLAGS}
            -DCMAKE_LINKER:FILEPATH=${CMAKE_LINKER}
            -DCMAKE_MAKE_PROGRAM:FILEPATH=${CMAKE_MAKE_PROGRAM}
            -DCMAKE_NM:FILEPATH=${CMAKE_NM}
            -DCMAKE_OBJCOPY:FILEPATH=${CMAKE_OBJCOPY}
            -DCMAKE_OBJDUMP:FILEPATH=${CMAKE_OBJDUMP}
            -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
            -DCMAKE_STRIP:FILEPATH=${CMAKE_STRIP})

    # On Windows, the necessary dlls get installed into the bin directory and
    # not lib.  This is the default behavior on CMake and I'm not sure how to
    # overwrite.  As such, we just forcefully move the files, so that the rest
    # of the scripts don't get messed up.
    if(WIN32 AND NOT EXISTS "${CMAKE_BINARY_DIR}/thirdparty/lib/${CMAKE_SHARED_LIBRARY_PREFIX}lapack${CMAKE_SHARED_LIBRARY_SUFFIX}")
        add_custom_target(
            fix_lapack_install_copy_bin
            ALL
            COMMAND ${CMAKE_COMMAND} -E copy_directory
                ${CMAKE_BINARY_DIR}/thirdparty/bin
                ${CMAKE_BINARY_DIR}/thirdparty/lib
            COMMAND ${CMAKE_COMMAND} -E remove_directory
                ${CMAKE_BINARY_DIR}/thirdparty/bin
            DEPENDS LAPACK)
    endif()

    # Set the location of BLAS and LAPACK
    set(BLAS_LIBRARY
        "${CMAKE_BINARY_DIR}/thirdparty/lib/${CMAKE_SHARED_LIBRARY_PREFIX}blas${CMAKE_SHARED_LIBRARY_SUFFIX}"
        CACHE FILEPATH "Path to a library." FORCE)
    set(BLAS_LIBRARIES ${BLAS_LIBRARY})
    mark_as_advanced(BLAS_LIBRARY)

    set(LAPACK_LIBRARY
        "${CMAKE_BINARY_DIR}/thirdparty/lib/${CMAKE_SHARED_LIBRARY_PREFIX}lapack${CMAKE_SHARED_LIBRARY_SUFFIX}"
        CACHE FILEPATH "Path to a library." FORCE)
    set(LAPACK_LIBRARIES ${LAPACK_LIBRARY})
    mark_as_advanced(LAPACK_LIBRARY)
else()
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
endif()

# Pass all the library variables to the parent
set(BLAS_LIBRARY ${BLAS_LIBRARY} PARENT_SCOPE)
set(BLAS_LIBRARIES ${BLAS_LIBRARIES} PARENT_SCOPE)

set(LAPACK_LIBRARY ${LAPACK_LIBRARY} PARENT_SCOPE)
set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} PARENT_SCOPE)

# Build jsonlab
if(ENABLE_MATLAB OR ENABLE_OCTAVE)
    mark_as_advanced(CLEAR ENABLE_BUILD_JSONLAB)
    set(ENABLE_BUILD_JSONLAB OFF CACHE BOOL "Build jsonlab from source?")
    if(ENABLE_BUILD_JSONLAB)
        # Determine the location of the jsonlab archive
        set(JSONLAB_ARCHIVE_UNDEFINED "JSONLAB_ARCHIVE-UNDEFINED")
        set(JSONLAB_ARCHIVE ${JSONLAB_ARCHIVE_UNDEFINED} CACHE FILEPATH
            "Jsonlab archive with full path.")
        # Make sure that it is defined and throw an error message if it is not
        if(JSONLAB_ARCHIVE STREQUAL JSONLAB_ARCHIVE_UNDEFINED)
            mark_as_advanced(CLEAR JSONLAB_ARCHIVE)
            message(FATAL_ERROR
                "The jsonlab archive with full path must be defined in the variable JSONLAB_ARCHIVE.")
        endif()
        mark_as_advanced(FORCE JSONLAB_ARCHIVE)

        include(ExternalProject)
        ExternalProject_Add(
            jsonlab
            URL ${JSONLAB_ARCHIVE}
            PATCH_COMMAND ${CMAKE_COMMAND} -E copy
                "${CMAKE_CURRENT_SOURCE_DIR}/jsonlab/CMakeLists.txt"
                CMakeLists.txt
            CMAKE_ARGS
                -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_BINARY_DIR}/thirdparty)

        # Set the location of jsonlab
        set(JSONLAB_DIR
            "${CMAKE_BINARY_DIR}/thirdparty/share/jsonlab"
            CACHE FILEPATH "Path to a library." FORCE)
        mark_as_advanced(JSONLAB_DIR)
    else()
        find_package(Jsonlab REQUIRED)
    endif()
endif()
