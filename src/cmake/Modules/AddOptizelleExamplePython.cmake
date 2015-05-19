# Installs a Python example 
macro(add_optizelle_example_python name)
   
    # Make sure we're installing examples
    if(ENABLE_PYTHON_EXAMPLES)
        # Install it to the standard location
        install(FILES "${name}.py"
            DESTINATION share/optizelle/examples/${name})
    endif()

endmacro()
