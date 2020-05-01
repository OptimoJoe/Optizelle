project(optizelle_python)

# Grab the current directory for some Python includes
set(OPTIZELLE_PYTHON_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})

# Compile the generic python optimization sources
set(optizelle_python_srcs "Utility.cpp")
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${OPTIZELLE_PYTHON_INCLUDE_DIRS})
add_library(optizelle_python MODULE ${optizelle_python_srcs})
target_link_libraries(optizelle_python optizelle_shared)
    
if(WIN32)
  target_link_libraries(optizelle_python optizelle_shared ${PYTHON_LIBRARIES})
elseif(APPLE)
  set_target_properties(optizelle_python PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif()

set_target_properties(optizelle_python PROPERTIES OUTPUT_NAME Utility)
set_target_properties(optizelle_python PROPERTIES PREFIX "")
if(APPLE)
    set_target_properties(optizelle_python PROPERTIES SUFFIX ".so")
elseif(WIN32)
    set_target_properties(optizelle_python PROPERTIES SUFFIX ".pyd")
endif()
if(APPLE)
    set_target_properties(optizelle_python PROPERTIES INSTALL_RPATH "@loader_path/../../../../lib;@loader_path/../../thirdparty/lib")
elseif(UNIX)
    set_target_properties(optizelle_python PROPERTIES INSTALL_RPATH "\$ORIGIN/../../../../lib;\$ORIGIN/../../thirdparty/lib")
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

if (WIN32)
  set (PYTHON_SITE_PACKAGES Lib/site-packages)
else ()
  set (PYTHON_SITE_PACKAGES lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages)
endif ()

# Compile the different problem classes and utility functions
add_subdirectory(Unconstrained)
add_subdirectory(InequalityConstrained)
add_subdirectory(EqualityConstrained)
add_subdirectory(Constrained)
add_subdirectory(json)

install(
    FILES __init__.py Exception.py Messaging.py Properties.py Functions.py
        Enumerated.py
    DESTINATION ${PYTHON_SITE_PACKAGES}/Optizelle)
install(
    TARGETS optizelle_python
    DESTINATION ${PYTHON_SITE_PACKAGES}/Optizelle)
