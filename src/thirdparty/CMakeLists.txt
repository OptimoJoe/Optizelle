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
        PATCH_COMMAND ${CMAKE_COMMAND} -E copy
            "${CMAKE_CURRENT_SOURCE_DIR}/jsoncpp/CMakeLists.txt"
            CMakeLists.txt
        CMAKE_ARGS
            -DCMAKE_BUILD_TYPE:STRING=RELEASE
            -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_BINARY_DIR}/thirdparty
    )

    # Set the location of jsoncpp
    set(JSONCPP_LIBRARY
        "${CMAKE_BINARY_DIR}/thirdparty/lib/${CMAKE_SHARED_LIBRARY_PREFIX}json${CMAKE_SHARED_LIBRARY_SUFFIX}"
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
    ExternalProject_Add(
        LAPACK 
        URL ${LAPACK_ARCHIVE}
        CMAKE_ARGS
            -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_BINARY_DIR}/thirdparty
            -DBUILD_TESTING:BOOL=OFF
            -DCMAKE_BUILD_TYPE:STRING=RELEASE
            -DBUILD_SHARED_LIBS:BOOL=ON
    )

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
    
# Install any packages that we compiled 
if(ENABLE_BUILD_JSONCPP OR ENABLE_BUILD_BLAS_AND_LAPACK)
    install(DIRECTORY ${CMAKE_BINARY_DIR}/thirdparty/ DESTINATION .)
endif()
