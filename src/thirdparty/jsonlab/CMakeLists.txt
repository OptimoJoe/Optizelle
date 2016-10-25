project(jsonlab)
cmake_minimum_required(VERSION 3.1.3)

# All we're doing is unpacking things and copying it over
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION share)
