project(peopt_cpp)

# Generate the proper headers for interfacing to BLAS and LAPACK
enable_language(Fortran)
include(FortranCInterface)
FortranCInterface_HEADER(FortranCInterface.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Turn on position independent code.  This simplifies  building both static
# and shared libraries.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Compile the rest of the library
set(peopt_cpp_srcs "vspaces.cpp" "peopt.cpp" "linalg.cpp")
add_library(peopt_cpp OBJECT ${peopt_cpp_srcs})
    
# Build the libraries 
add_library(peopt_static STATIC $<TARGET_OBJECTS:peopt_cpp>)
add_library(peopt_shared SHARED $<TARGET_OBJECTS:peopt_cpp>)
add_library(peopt_shared_complete SHARED $<TARGET_OBJECTS:peopt_cpp>)
target_link_libraries(peopt_shared_complete
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

# This is required for Windows to see an import library, which is
# required to build a dll.
target_link_libraries(peopt_shared peopt_shared_complete)

# Install the headers associated with the library.  
install(FILES
    vspaces.h
    peopt.h
    json.h
    linalg.h
    DESTINATION include/peopt)
install(TARGETS
    peopt_static
    peopt_shared
    peopt_shared_complete
    DESTINATION lib)    
