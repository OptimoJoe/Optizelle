# Sets up and runs Optizelle unit tests 
macro(add_optizelle_test_python name)

    # Make sure that tests are enabled
    if(ENABLE_PYTHON_UNIT)
        add_test("Execution_of_python_${name}" ${PYTHON_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/${name}.py")
        set_tests_properties("Execution_of_python_${name}"
            PROPERTIES ENVIRONMENT
                "PYTHONPATH=${CMAKE_BINARY_DIR}/src/python")
    endif()

endmacro()
