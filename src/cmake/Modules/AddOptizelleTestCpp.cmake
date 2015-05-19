# Sets up and runs Optizelle unit tests 
macro(add_optizelle_test_cpp name)

    # Make sure that tests are enabled
    if(ENABLE_CPP_UNIT)
        add_test("Execution_of_cpp_${name}" ${name})
    endif()

endmacro()
