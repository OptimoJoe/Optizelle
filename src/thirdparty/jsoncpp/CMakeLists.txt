project(jsoncpp)
cmake_minimum_required(VERSION 3.1.3)
enable_language(CXX)

# Set the source file prefix
set(source_prefix ${CMAKE_CURRENT_SOURCE_DIR}/src/lib_json/)

# Save the include directory
set(JSONCPP_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Grab all the sources for the library
set(jsoncpp_srcs
    "${source_prefix}/json_reader.cpp"
    "${source_prefix}/json_value.cpp"
    "${source_prefix}/json_writer.cpp"
)

# Locate the headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Compile everything
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
add_library(jsoncpp_object OBJECT ${jsoncpp_srcs})
add_library(jsoncpp_static STATIC $<TARGET_OBJECTS:jsoncpp_object> )
add_library(jsoncpp_shared SHARED $<TARGET_OBJECTS:jsoncpp_object> )
set_target_properties(jsoncpp_shared jsoncpp_static
    PROPERTIES OUTPUT_NAME json)

# Set the compiler to C++14
set_property(TARGET jsoncpp_object PROPERTY CXX_STANDARD 14)
set_property(TARGET jsoncpp_object PROPERTY CXX_STANDARD_REQUIRED ON)
set_property(TARGET jsoncpp_static PROPERTY CXX_STANDARD 14)
set_property(TARGET jsoncpp_static PROPERTY CXX_STANDARD_REQUIRED ON)
set_property(TARGET jsoncpp_shared PROPERTY CXX_STANDARD 14)
set_property(TARGET jsoncpp_shared PROPERTY CXX_STANDARD_REQUIRED ON)

# Install the libraries and headers
install(TARGETS jsoncpp_static jsoncpp_shared DESTINATION lib)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/json DESTINATION include)
