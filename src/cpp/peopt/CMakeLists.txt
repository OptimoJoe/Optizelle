cmake_minimum_required(VERSION 2.8.9)
project(peopt_cpp)

set(peopt_cpp_srcs "linalg.cpp")
include(FortranCInterface)
FortranCInterface_HEADER(${CMAKE_CURRENT_BINARY_DIR}/FortranCInterface.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})
add_library(peopt_cpp OBJECT ${peopt_cpp_srcs})
# Install the headers associated with the library.  In theory, json.h and
# peopt.h should already be installed
install(FILES linalg.h vspaces.h DESTINATION include/peopt)
