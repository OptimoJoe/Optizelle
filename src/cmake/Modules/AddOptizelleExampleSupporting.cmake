# Installs the files required for an Optizelle example 
macro(add_optizelle_example_supporting name)

    # Make sure that examples are enabled 
    if(ENABLE_CPP_EXAMPLES OR ENABLE_PYTHON_EXAMPLES OR ENABLE_MATLAB_EXAMPLES)

        # Grab the file list 
        set(files "${ARGN}")

        # Install the example files to the correct location
        install(FILES ${files}
            DESTINATION share/optizelle/examples/${name})
    endif()

endmacro()

macro(add_optizelle_example_supporting_restart name)

    # Make sure that examples are enabled 
    if(ENABLE_CPP_EXAMPLES OR ENABLE_PYTHON_EXAMPLES OR ENABLE_MATLAB_EXAMPLES)

        # Grab the file list 
        set(files "${ARGN}")

        # Install the example restart files to the correct location
        install(FILES ${files}
            DESTINATION share/optizelle/examples/${name}/restart)
    endif()

    # If we're doing unit tests, make sure the restart files are in the right
    # location 
    if(ENABLE_CPP_UNIT OR ENABLE_PYTHON_UNIT OR ENABLE_MATLAB_UNIT)
        add_custom_target(
            ${name}_restart_files
            ALL
            COMMAND ${CMAKE_COMMAND} -E copy_directory
                ${CMAKE_SOURCE_DIR}/src/examples/${name}/restart
                ${CMAKE_BINARY_DIR}/src/examples/${name}/restart)
endif()

endmacro()
