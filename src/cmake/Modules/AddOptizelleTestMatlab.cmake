# Sets up and runs Optizelle unit tests 
macro(add_optizelle_test_matlab name)

    # Make sure that tests are enabled
    if(ENABLE_MATLAB_UNIT)
        add_test( "Execution_of_matlab_${name}"
            ${MATLAB_EXECUTABLE} ${MATLAB_RUN_FLAG} 
            "${executable}${name}()")
        set_tests_properties("Execution_of_matlab_${name}"
            PROPERTIES ENVIRONMENT
                "MATLABPATH=${CMAKE_BINARY_DIR}/src/matlab:${CMAKE_CURRENT_SOURCE_DIR}")
        set_tests_properties("Execution_of_matlab_${name}"
            PROPERTIES ENVIRONMENT
                "OCTAVE_PATH=${CMAKE_BINARY_DIR}/src/matlab:${CMAKE_CURRENT_SOURCE_DIR}")
    endif()

endmacro()
