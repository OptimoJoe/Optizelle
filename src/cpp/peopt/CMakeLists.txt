project(peopt_cpp)

# Generate the proper headers for interfacing to BLAS and LAPACK
include(FortranCInterface)
FortranCInterface_HEADER(FortranCInterface.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Compile the rest of the library
set(peopt_cpp_srcs "linalg.cpp")
add_library(peopt_cpp OBJECT ${peopt_cpp_srcs})

# Install the headers associated with the library.  In theory, json.h and
# peopt.h should already be installed
install(FILES vspaces.h DESTINATION include/peopt)
