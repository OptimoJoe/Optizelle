# Enable C++
enable_language(CXX)

# Figure out if we should look for OpenMP
mark_as_advanced(CLEAR ENABLE_OPENMP)
set(ENABLE_OPENMP OFF CACHE BOOL "Enable OpenMP?")
if(ENABLE_OPENMP)
    # If we enable OpenMP, go ahead and make it required and find the
    # library.
    find_package(OpenMP REQUIRED)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" PARENT_SCOPE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" PARENT_SCOPE)
endif()

# Set the Optizelle include directories
set(OPTIZELLE_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})
set(OPTIZELLE_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)

# Compile the Optizelle core 
add_subdirectory(optizelle)
