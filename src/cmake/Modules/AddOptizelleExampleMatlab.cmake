# Installs a Matlab example 
macro(add_optizelle_example_matlab name)
   
    # Make sure we're installing examples
    if(ENABLE_MATLAB_EXAMPLES)
        # Install it to the standard location
        install(FILES "${name}.m"
            DESTINATION share/optizelle/examples/${name})
    endif()

endmacro()
