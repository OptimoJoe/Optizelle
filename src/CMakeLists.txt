# Add a module directory in case we need to find some packages
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
    "${CMAKE_SOURCE_DIR}/src/cmake/Modules/")
    
# Add some useful macros for examples and unit testing 
include(Optizelle)

# Grab our third party library dependencies
add_subdirectory(thirdparty)

# Compile our C++ pieces 
add_subdirectory(cpp)

# Compile our Python pieces
add_subdirectory(python)

# Compile our MATLAB/Octave pieces
add_subdirectory(matlab_octave)

# Compile our unit tests
add_subdirectory(unit)

# Compile our examples 
add_subdirectory(examples)

# See if we're going to build the documentation for the package
set(ENABLE_DOCUMENTATION OFF CACHE BOOL "Enable documentation generation?")
if(ENABLE_DOCUMENTATION)
    add_subdirectory(doc)

# If we're not building the documentation, hide some of the parameters
else()
    set(ENABLE_A4_PAPER OFF CACHE BOOL "Generate A4 sized documentation?  If false, the document is letter sized." FORCE)
    mark_as_advanced(FORCE ENABLE_A4_PAPER)
endif()
