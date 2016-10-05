# Check if we want to build the MATLAB interface 
mark_as_advanced(CLEAR ENABLE_MATLAB)
set(ENABLE_MATLAB OFF CACHE BOOL "Enable the MATLAB build?")

# Check if we want to build the Octave interface 
mark_as_advanced(CLEAR ENABLE_OCTAVE)
set(ENABLE_OCTAVE OFF CACHE BOOL "Enable the Octave build?")

# MATLAB specific options
if(ENABLE_MATLAB)
    # Find the MATLAB libraries 
    find_package(MATLAB REQUIRED)

    # Install the setup file
    install(FILES setupOptizelle.m DESTINATION share/optizelle/matlab/)

    # If we're on Windows, install a readme file that explains how to set the
    # MATLAB path
    if(WIN32)
        install(FILES readme_windows.txt RENAME readme.txt DESTINATION share/optizelle/matlab/)
    endif()
endif()

# Octave specific options
if(ENABLE_OCTAVE)
    # Find the Octave libraries 
    find_package(Octave REQUIRED)

    # Install the setup file
    install(FILES setupOptizelle.m DESTINATION share/optizelle/octave/)
endif()

# Options for both
if(ENABLE_MATLAB OR ENABLE_OCTAVE)
    # Common includes
    include_directories(${OPTIZELLE_INCLUDE_DIRS})
    include_directories(${JSONCPP_INCLUDE_DIRS})

    # Add the directory with the sources
    add_subdirectory(optizelle)
endif()
