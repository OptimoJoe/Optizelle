cmake_minimum_required(VERSION 2.8.9)
project(peopt_python)

set(peopt_python_srcs "peopt.cpp")
add_library(peopt_python OBJECT ${peopt_python_srcs})
