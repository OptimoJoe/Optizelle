# Sets up and runs tests on sparse SDPA formatted files.
macro(add_optizelle_sdpa_cpp dats params executable)

    # Make sure that tests are enabled
    if(ENABLE_CPP_UNIT)

        # Grab all of the .json paramters files and the .dat-s files
        file(GLOB_RECURSE dats ${CMAKE_CURRENT_SOURCE_DIR} ${dats})
        list(SORT dats)
        file(GLOB_RECURSE params ${CMAKE_CURRENT_SOURCE_DIR} ${params})
        list(SORT params)


        # We start indexing from 0, but the length parameter starts at1.
        # We do the math here to fix it.
        list(LENGTH dats len1)
        math(EXPR len2 "${len1} - 1")

        # Loop over each of the tests
        foreach(index RANGE ${len2})
            # Grab the particular unit test
            list(GET dats ${index} dat)

            # Separate out the phase-1 and phase-2 parameters
            math(EXPR phase1_idx "2 * ${index}")
            math(EXPR phase2_idx "2 * ${index} + 1")
            list(GET params ${phase1_idx} phase1)
            list(GET params ${phase2_idx} phase2)

            # Run the optimization
            add_test("Execution_of_${dat}"
                ${executable} ${dat} ${phase1} ${phase2})
            
            # Diff the result of the optimization against the known solution
            add_test("Solution_to_${dat}"
                "${CMAKE_BINARY_DIR}/src/unit/utility/diff_restart"
                ${phase2} solution.json)
        endforeach()
    endif()
endmacro()
