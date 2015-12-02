# Sets up and runs Optizelle unit tests 
macro(add_optizelle_test_cpp name)

    # Make sure that tests are enabled
    if(ENABLE_CPP_UNIT)
        # Grab the absolute path of the test
        get_filename_component(fullname ${name} ABSOLUTE)

        # Add the test
        add_test("Execution_of_cpp_${fullname}" ${name})
    endif()

endmacro()
