# Sets up and runs Optizelle unit tests 
macro(add_optizelle_test_cpp files executable)

    # Make sure that tests are enabled
    if(ENABLE_CPP_UNIT)

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
            add_test("Execution_of_${unit}" ${executable} ${unit})

            # Diff the result of the optimization against the known solution
            add_test("Solution_to_${unit}"
                "${CMAKE_BINARY_DIR}/src/cpp/unit/utility/diff_restart"
                ${unit} solution.json)
        endforeach()
    endif()

endmacro()
