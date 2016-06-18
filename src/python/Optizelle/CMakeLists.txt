project(optizelle_python)

# Grab the current directory for some Python includes
set(OPTIZELLE_PYTHON_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})

# Compile the generic python optimization sources
set(optizelle_python_srcs "Utility.cpp")
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${OPTIZELLE_PYTHON_INCLUDE_DIRS})
add_library(optizelle_python SHARED ${optizelle_python_srcs})
target_link_libraries(optizelle_python
    optizelle_shared
    ${PYTHON_LIBRARIES})
set_target_properties(optizelle_python PROPERTIES OUTPUT_NAME Utility)
set_target_properties(optizelle_python PROPERTIES PREFIX "")
if(APPLE)
    set_target_properties(optizelle_python PROPERTIES SUFFIX ".so") 
elseif(WIN32)
    set_target_properties(optizelle_python PROPERTIES SUFFIX ".pyd")
endif()

# Copy in the rest of the Python files for the unit tests.  Otherwise, we
# end of dirtying our source tree.
if(ENABLE_PYTHON_UNIT)
    add_custom_target(
        python_files
        ALL
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${CMAKE_SOURCE_DIR}/src/python/Optizelle
            ${CMAKE_BINARY_DIR}/src/python/Optizelle)
endif()

# Compile the different problem classes and utility functions
add_subdirectory(Unconstrained)
add_subdirectory(InequalityConstrained)
add_subdirectory(EqualityConstrained)
add_subdirectory(Constrained)
add_subdirectory(json)

install(
    FILES __init__.py Exception.py Messaging.py Properties.py Functions.py
        Enumerated.py
    DESTINATION share/optizelle/python/Optizelle)
install(
    TARGETS optizelle_python 
    DESTINATION share/optizelle/python/Optizelle)
