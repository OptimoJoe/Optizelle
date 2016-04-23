# Runs a MATLAB/Octave json test 
macro(run_matlab_octave_json_test files script type executable flags environment validated)
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
        add_test(
            "Execution_of_${type}_${unit}"
            ${executable}
            ${flags}
            "${script}('${unit}'),exit")

        # Make sure that Optizelle is available 
        set_tests_properties(
            "Execution_of_${type}_${unit}"
            PROPERTIES ENVIRONMENT
            ${environment})

        # Diff the result of the optimization against the known solution
        if(validated)
            add_test("Solution_to_${type}_${unit}"
                "${CMAKE_BINARY_DIR}/src/unit/utility/diff_restart"
                ${unit} solution.json)
        endif()
    endforeach()
endmacro()

# Sets up and runs Optizelle unit tests 
macro(add_optizelle_json_test_matlab files script validated)
    # MATLAB on Windows does not set the path through an environment variable,
    # so we can't easily run the unit tests
    if(ENABLE_MATLAB_UNIT)
        run_matlab_octave_json_test(
            ${files}
            ${script}
            matlab
            ${MATLAB_EXECUTABLE}
            "-nosplash -nodesktop -r"
            "MATLABPATH=${CMAKE_BINARY_DIR}/src/matlab:${CMAKE_CURRENT_SOURCE_DIR}"
            validated)
    endif()

    # Run the Octave tests
    if(ENABLE_OCTAVE_UNIT)
        run_matlab_octave_json_test(
            ${files}
            ${script}
            octave 
            ${OCTAVE_EXECUTABLE}
            "--eval"
            "OCTAVE_PATH=${CMAKE_BINARY_DIR}/src/octave:${CMAKE_CURRENT_SOURCE_DIR}"
            validated)
    endif()
endmacro()
