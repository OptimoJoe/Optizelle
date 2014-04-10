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
