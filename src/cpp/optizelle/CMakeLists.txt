project(optizelle_cpp)

# Generate the proper headers for interfacing to BLAS and LAPACK
enable_language(Fortran)
include(FortranCInterface)
FortranCInterface_HEADER(FortranCInterface.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Grab the rest of the include files
include_directories(${OPTIZELLE_INCLUDE_DIRS})
include_directories(${JSONCPP_INCLUDE_DIRS})

# Utility sources.  We use these routines for things like preprocessing and
# other pieces that don't need all of Optizelle, so compile them separately.
set(utility_srcs "exception.cpp" "stream.cpp")
add_library(utility OBJECT ${utility_srcs})

# Compile the core of Optizelle
set(optizelle_cpp_srcs "vspaces.cpp" "optizelle.cpp" "linalg.cpp" "json.cpp")
add_library(optizelle_cpp OBJECT ${optizelle_cpp_srcs})

# Package everything together
add_library(optizelle_static STATIC
    $<TARGET_OBJECTS:utility>
    $<TARGET_OBJECTS:optizelle_cpp>)
set_target_properties(optizelle_static PROPERTIES OUTPUT_NAME optizelle)

add_library(optizelle_shared SHARED
    $<TARGET_OBJECTS:utility>
    $<TARGET_OBJECTS:optizelle_cpp>)
set_target_properties(optizelle_shared PROPERTIES OUTPUT_NAME optizelle)
target_link_libraries(optizelle_shared
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})
if(APPLE)
    set_target_properties(optizelle_shared PROPERTIES INSTALL_RPATH "@loader_path/../share/optizelle/thirdparty/lib")
elseif(UNIX)
    set_target_properties(optizelle_shared PROPERTIES INSTALL_RPATH "\$ORIGIN/../share/optizelle/thirdparty/lib")
endif()

# Install the headers associated with the library.
install(FILES
    vspaces.h
    optizelle.h
    json.h
    linalg.h
    exception.h
    stream.h
    DESTINATION include/optizelle)
install(TARGETS
    optizelle_static
    optizelle_shared
    DESTINATION lib)

# Grab required runtime
if(MINGW)
    get_filename_component(Mingw_Path ${CMAKE_CXX_COMPILER} PATH)
    set(CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS
        ${Mingw_Path}/libwinpthread-1.dll
        ${Mingw_Path}/libgcc_s_seh-1.dll
        ${Mingw_Path}/libstdc++-6.dll
        ${Mingw_Path}/libgfortran-4.dll
        ${Mingw_Path}/libquadmath-0.dll
        ${Mingw_Path}/libgomp-1.dll)
    file(
	COPY ${CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS}
	DESTINATION ${CMAKE_BINARY_DIR}/thirdparty/lib)
elseif(APPLE AND CPACK_DRAGNDROP)
    get_filename_component(Compiler_Path ${CMAKE_Fortran_COMPILER} PATH)
    set(Lib_Path "${Compiler_Path}/../lib/libgcc")
    file(COPY
        ${Lib_Path}/libgfortran.3.dylib
        ${Lib_Path}/libquadmath.0.dylib
        ${Lib_Path}/libgcc_s.1.dylib
        DESTINATION
        ${CMAKE_BINARY_DIR}/thirdparty/lib)
endif()
