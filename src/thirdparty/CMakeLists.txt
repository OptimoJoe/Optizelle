project(third_party_libraries)

# Build jsoncpp
mark_as_advanced(CLEAR ENABLE_BUILD_JSONCPP)
set(ENABLE_BUILD_JSONCPP OFF CACHE BOOL "Build jsoncpp from source?")
if(ENABLE_BUILD_JSONCPP)
    include(ExternalProject)
    ExternalProject_Add(
        JsonCpp
        URL ${CMAKE_BINARY_DIR}/jsoncpp-src-0.6.0-rc2.tar.gz
        URL_MD5 363e2f4cbd3aeb63bf4e571f377400fb
        PATCH_COMMAND ${CMAKE_COMMAND} -E copy
            "${CMAKE_CURRENT_SOURCE_DIR}/jsoncpp/CMakeLists.txt"
            CMakeLists.txt
        CMAKE_ARGS
            -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_BINARY_DIR}/installed
            -DCMAKE_BUILD_TYPE:STRING=RELEASE
    )

    # Set the location of jsoncpp
    set(JSONCPP_LIBRARY
        "${CMAKE_BINARY_DIR}/installed/lib/${CMAKE_STATIC_LIBRARY_PREFIX}json${CMAKE_STATIC_LIBRARY_SUFFIX}"
        CACHE FILEPATH "Path to a library." FORCE)
    set(JSONCPP_LIBRARIES ${JSONCPP_LIBRARY}) 
    set(JSONCPP_INCLUDE_DIR ${CMAKE_BINARY_DIR}/installed/include CACHE
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
    include(ExternalProject)
    ExternalProject_Add(
        LAPACK 
        URL ${CMAKE_BINARY_DIR}/lapack-3.4.2.tgz
        URL_MD5 61bf1a8a4469d4bdb7604f5897179478 
        CMAKE_ARGS
            -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_BINARY_DIR}/installed
            -DBUILD_TESTING:BOOL=OFF
            -DCMAKE_BUILD_TYPE:STRING=RELEASE
    )

    # Set the location of BLAS and LAPACK 
    set(BLAS_LIBRARY
        "${CMAKE_BINARY_DIR}/installed/lib/${CMAKE_STATIC_LIBRARY_PREFIX}blas${CMAKE_STATIC_LIBRARY_SUFFIX}"
        CACHE FILEPATH "Path to a library." FORCE)
    set(BLAS_LIBRARIES ${BLAS_LIBRARY}) 
    mark_as_advanced(BLAS_LIBRARY)
    
    set(LAPACK_LIBRARY
        "${CMAKE_BINARY_DIR}/installed/lib/${CMAKE_STATIC_LIBRARY_PREFIX}lapack${CMAKE_STATIC_LIBRARY_SUFFIX}"
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
