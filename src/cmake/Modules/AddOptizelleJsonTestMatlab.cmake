# Sets up and runs Optizelle unit tests 
macro(add_optizelle_json_test_matlab files executable)
    
    # Grab any extra arguments
    set (extra_macro_args ${ARGN})
    list(LENGTH extra_macro_args num_extra_args)
    if (${num_extra_args} GREATER 0)
        list(GET extra_macro_args 0 optional_arg)
    endif()

    # Make sure that tests are enabled
    if(ENABLE_MATLAB_UNIT)

        # Grab all of the different .json test files
        file(GLOB_RECURSE units ${CMAKE_CURRENT_SOURCE_DIR} ${files})

        # Make sure we run everything in order
        list(SORT units)

        # Remove any optional json file from this list
        if (${num_extra_args} GREATER 0)
            list(REMOVE_ITEM units ${optional_arg})
        endif()

        # We start indexing from 0, but the length parameter starts at1.
        # We do the math here to fix it.
        list(LENGTH units len1)
        math(EXPR len2 "${len1} - 1")

        # Loop over each of the tests
        foreach(index RANGE ${len2})

            # Grab the particular unit test
            list(GET units ${index} unit)
                
            # Figure out whether we have MATLAB or Octave.  If is_matlab is
            # negative we have Octave, otherwise we have MATLAB. 
            string(TOLOWER ${MATLAB_EXECUTABLE} prog)
            string(FIND ${prog} "matlab" is_matlab)

            # Run the optimization
            if (${num_extra_args} EQUAL 0)
                if(${is_matlab} LESS 0)
                    add_test( "Execution_of_matlab_${unit}"
                        ${MATLAB_EXECUTABLE} "--eval" 
                        "${executable}('${unit}'),exit")
                else()
                    add_test( "Execution_of_matlab_${unit}"
                        ${MATLAB_EXECUTABLE} "-nosplash -nodesktop -r"
                        "\"${executable}('${unit}'),exit\"")
                endif()
            else()
                list(GET extra_macro_args 0 optional_arg)
                if(${is_matlab} LESS 0)
                    add_test( "Execution_of_matlab_${unit}"
                        ${MATLAB_EXECUTABLE} "--eval" 
                        "${executable}('${unit}','${optional_arg}'),exit")
                else()
                    add_test( "Execution_of_matlab_${unit}"
                        ${MATLAB_EXECUTABLE} "-nosplash -nodesktop -r"
                        "\"${executable}('${unit}','${optional_arg}'),exit\"")
                endif()
            endif()

            set_tests_properties("Execution_of_matlab_${unit}"
                PROPERTIES ENVIRONMENT
                    "MATLABPATH=${CMAKE_BINARY_DIR}/src/matlab:${CMAKE_CURRENT_SOURCE_DIR}")
            set_property(TEST "Execution_of_matlab_${unit}"
                APPEND PROPERTY ENVIRONMENT
                    "OCTAVE_PATH=${CMAKE_BINARY_DIR}/src/matlab:${CMAKE_CURRENT_SOURCE_DIR}")

            # Diff the result of the optimization against the known solution
            add_test("Solution_to_matlab_${unit}"
                "${CMAKE_BINARY_DIR}/src/unit/utility/diff_restart"
                ${unit} solution.json)
        endforeach()
    endif()

endmacro()
