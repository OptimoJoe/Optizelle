# Sets up and runs Optizelle unit tests 
macro(add_optizelle_json_test_matlab files executable)

    # Make sure that tests are enabled
    if(ENABLE_MATLAB_UNIT)

        # Grab all of the different .json test files
        file(GLOB_RECURSE units ${CMAKE_CURRENT_SOURCE_DIR} ${files})

        # Make sure we run everything in order
        list(SORT units)

        # We start indexing from 0, but the length parameter starts at1.
        # We do the math here to fix it.
        list(LENGTH units len1)
        math(EXPR len2 "${len1} - 1")

        # Loop over each of the tests
        foreach(index RANGE ${len2})

            # Grab the particular unit test
            list(GET units ${index} unit)

            # Run the optimization
            add_test( "Execution_of_matlab_${unit}"
                ${MATLAB_EXECUTABLE} ${MATLAB_RUN_FLAG} 
                "${executable}('${unit}')")
            set_tests_properties("Execution_of_matlab_${unit}"
                PROPERTIES ENVIRONMENT
                    "MATLABPATH=${CMAKE_BINARY_DIR}/src/matlab:${CMAKE_CURRENT_SOURCE_DIR}")
            set_tests_properties("Execution_of_matlab_${unit}"
                PROPERTIES ENVIRONMENT
                    "OCTAVE_PATH=${CMAKE_BINARY_DIR}/src/matlab:${CMAKE_CURRENT_SOURCE_DIR}")

            # Diff the result of the optimization against the known solution
            add_test("Solution_to_matlab_${unit}"
                "${CMAKE_BINARY_DIR}/src/unit/utility/diff_restart"
                ${unit} solution.json)
        endforeach()
    endif()

endmacro()
