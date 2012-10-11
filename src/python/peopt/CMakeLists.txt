cmake_minimum_required(VERSION 2.8.9)
project(peopt_python)

set(peopt_python_srcs "peopt.cpp")
include_directories(${PYTHON_INCLUDE_DIRS} ${JSONCPP_INCLUDE_DIRS})
add_library(peopt_python OBJECT ${peopt_python_srcs})
install(FILES peopt.py DESTINATION share/peopt/python)
