project(optizelle_cpp)

# Generate the proper headers for interfacing to BLAS and LAPACK
enable_language(Fortran)
include(FortranCInterface)
FortranCInterface_HEADER(FortranCInterface.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Grab the rest of the include files
include_directories(${OPTIZELLE_INCLUDE_DIRS})
include_directories(${JSONCPP_INCLUDE_DIRS})

# Compile the library
set(optizelle_cpp_srcs "vspaces.cpp" "optizelle.cpp" "linalg.cpp" "json.cpp")
add_library(optizelle_cpp OBJECT ${optizelle_cpp_srcs})
    
# Package everything together 
add_library(optizelle_static STATIC $<TARGET_OBJECTS:optizelle_cpp>)
set_target_properties(optizelle_static PROPERTIES OUTPUT_NAME optizelle)

add_library(optizelle_shared SHARED $<TARGET_OBJECTS:optizelle_cpp>)
set_target_properties(optizelle_shared PROPERTIES OUTPUT_NAME optizelle)
target_link_libraries(optizelle_shared
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

# Install the headers associated with the library.  
install(FILES
    vspaces.h
    optizelle.h
    json.h
    linalg.h
    DESTINATION include/optizelle)
install(TARGETS
    optizelle_static
    optizelle_shared
    DESTINATION lib)    
