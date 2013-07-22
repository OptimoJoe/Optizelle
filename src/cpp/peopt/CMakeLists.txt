project(peopt_cpp)

# Generate the proper headers for interfacing to BLAS and LAPACK
include(FortranCInterface)
FortranCInterface_HEADER(FortranCInterface.h)

# Grabs the peopt headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)

# Grabs the FortranCInterface header
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Compile the rest of the library
set(peopt_cpp_srcs "peopt.cpp" "linalg.cpp")
add_library(peopt_cpp OBJECT ${peopt_cpp_srcs})

# Install the headers associated with the library.  
install(FILES
    vspaces.h
    peopt.h
    json.h
    linalg.h
    DESTINATION include/peopt)
