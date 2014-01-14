project(optizelle_python)

# Turn on position independent code.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Grab the current directory for some Python includes
set(OPTIZELLE_PYTHON_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})

# Compile the generic python optimization sources
set(optizelle_python_srcs "Utility.cpp")
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${OPTIZELLE_PYTHON_INCLUDE_DIRS})
add_library(optizelle_python SHARED ${optizelle_python_srcs})
target_link_libraries(optizelle_python
    optizelle_static
    ${PYTHON_LIBRARIES}
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})
set_target_properties(optizelle_python
    PROPERTIES OUTPUT_NAME Utility)
set_target_properties(optizelle_python
    PROPERTIES PREFIX "") 

# Copy in the rest of the Python files for the unit tests.  Otherwise, we
# end of dirtying our source tree.
add_custom_command(
    TARGET optizelle_python PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/src/python/Optizelle
        $<TARGET_FILE_DIR:optizelle_python>)

# Compile the different problem classes and utility functions
add_subdirectory(Unconstrained)
add_subdirectory(InequalityConstrained)
add_subdirectory(EqualityConstrained)
add_subdirectory(Constrained)
add_subdirectory(json)

install(
    FILES __init__.py
    DESTINATION share/optizelle/python/Optizelle)
install(
    TARGETS optizelle_python 
    DESTINATION share/optizelle/python/Optizelle)
