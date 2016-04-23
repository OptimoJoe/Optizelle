# Adds an example both for the unit test as well as an exmaple
macro(add_example name interfaces supporting)
    # Set the base installation directory
    set(basedir "share/optizelle/examples/")

    # Grab and install the problem setups
    file(GLOB_RECURSE setups ${CMAKE_CURRENT_SOURCE_DIR} "*.json")
    install(FILES ${setups} DESTINATION share/optizelle/examples/${name})

    # Installs any additional supporting files
    install(FILES ${supporting} DESTINATION share/optizelle/examples/${name})

    # Compile and install the particular problem interface 
    foreach(interface ${interfaces})
        if(${interface} STREQUAL "cpp")

            # Compile and link
            if(ENABLE_CPP_EXAMPLES OR ENABLE_CPP_UNIT)
                include_directories(${OPTIZELLE_INCLUDE_DIRS})
                include_directories(${JSONCPP_INCLUDE_DIRS})
                add_executable(${name} "${name}.cpp")
                target_link_libraries(${name}
                    optizelle_static
                    ${JSONCPP_LIBRARIES}
                    ${LAPACK_LIBRARIES}
                    ${BLAS_LIBRARIES})
            endif()
                
            # Installation
            if(ENABLE_CPP_EXAMPLES)
                install(FILES "${name}.cpp" DESTINATION "${basedir}/${name}")
                install(TARGETS ${name} DESTINATION "${basedir}/${name}")
            endif()

        elseif(${interface} STREQUAL "python")
            # Installation
            if(ENABLE_PYTHON_EXAMPLES)
                install(FILES "${name}.py" DESTINATION "${basedir}/${name}")
            endif()

        elseif(${interface} STREQUAL "matlab" OR ${interface} STREQUAL "octave")
            # Installation 
            if(ENABLE_MATLAB_EXAMPLES OR ENABLE_OCTAVE_EXAMPLES)
                install(FILES "${name}.m" DESTINATION "${basedir}/${name}")
            endif()
        endif()
    endforeach()
endmacro()

# Adds a unit test, which may or may not validate the results
macro(add_unit name interfaces units validated)
    # Figure out the number of unit tests
    list(LENGTH units nunits)

    # If we have an extra argument, we grab it, at least the first one
    set(extra ${ARGN})
    list(LENGTH extra nextra)
    if(${nextra} GREATER 0)
        list(GET extra 0 extra)
    endif()

    # When we specify a unit test based on a json file
    if(nunits GREATER 0)
        # Loop over each of the tests
        foreach(unit ${units})
            # Grab the unit name
            get_filename_component(uname ${unit} NAME_WE)

            # Loop over each of the interface interfaces
            foreach(interface ${interfaces})
                if(${interface} STREQUAL "cpp" AND ENABLE_CPP_UNIT)
                    # Run the test
                    add_test("Execution_of_cpp_${name}_${uname}"
                        ${name}
                        ${unit}
                        ${extra})

                elseif(${interface} STREQUAL "python" AND ENABLE_PYTHON_UNIT)
                    # Run the test
                    add_test("Execution_of_python_${name}_${uname}"
                        ${PYTHON_EXECUTABLE}
                        "${CMAKE_CURRENT_SOURCE_DIR}/${name}.py"
                        ${unit}
                        ${extra})

                    # Make sure Optizelle is available
                    set_tests_properties("Execution_of_python_${name}_${uname}"
                        PROPERTIES ENVIRONMENT
                            "PYTHONPATH=${CMAKE_BINARY_DIR}/src/python")

                # Run the MATLAB test.  Note, this doesn't work on Windows
                # since we can't set an environment variable to change the
                # path.
                elseif(${interface} STREQUAL "matlab"
                    AND ENABLE_MATLAB_UNIT
                    AND NOT WIN32
                )
                    # If we have an extra argument, append it
                    if(${nextra} GREATER 0)
                        set(matlab_unit "${unit}','${extra}")
                    else()
                        set(matlab_unit "${unit}")
                    endif()

                    # Run the test
                    add_test(
                        "Execution_of_matlab_${name}_${uname}"
                        ${MATLAB_EXECUTABLE}
                        "-nosplash -nodesktop -r"
                        "\"${name}('${matlab_unit}'),exit\"")

                    # Make sure Optizelle is available 
                    set_tests_properties(
                        "Execution_of_matlab_${name}_${uname}"
                        PROPERTIES ENVIRONMENT
                        "MATLABPATH=${CMAKE_BINARY_DIR}/src/matlab:${CMAKE_CURRENT_SOURCE_DIR}")
                elseif(${interface} STREQUAL "octave" AND ENABLE_OCTAVE_UNIT)
                    # If we have an extra argument, append it
                    if(${nextra} GREATER 0)
                        set(octave_unit "${unit}','${extra}")
                    else()
                        set(octave_unit "${unit}")
                    endif()

                    # Run the test
                    add_test(
                        "Execution_of_octave_${name}_${uname}"
                        ${OCTAVE_EXECUTABLE}
                        "--eval"
                        "${name}('${octave_unit}'),exit")

                    # Make sure Optizelle is available 
                    set_tests_properties(
                        "Execution_of_octave_${name}_${uname}"
                        PROPERTIES ENVIRONMENT
                        "OCTAVE_PATH=${CMAKE_BINARY_DIR}/src/octave:${CMAKE_CURRENT_SOURCE_DIR}")
                endif()
                
                # Diff the result of the optimization against the known solution
                if(${validated})
                    add_test("Solution_to_${interface}_${name}_${uname}"
                        "${CMAKE_BINARY_DIR}/src/unit/utility/diff_restart"
                        ${unit} solution.json)
                endif()
            endforeach()
        endforeach()

    # When we're just running a script with no arguments
    else()
        # Loop over each of the interface interfaces
        foreach(interface ${interfaces})
            if(${interface} STREQUAL "cpp" AND ENABLE_CPP_UNIT)
                # Run the test
                add_test("Execution_of_cpp_${name}" ${name})

            elseif(${interface} STREQUAL "python" AND ENABLE_PYTHON_UNIT)
                # Run the test
                add_test("Execution_of_python_${name}" ${PYTHON_EXECUTABLE}
                    "${CMAKE_CURRENT_SOURCE_DIR}/${name}.py")

                # Make sure Optizelle is available
                set_tests_properties("Execution_of_python_${name}"
                    PROPERTIES ENVIRONMENT
                        "PYTHONPATH=${CMAKE_BINARY_DIR}/src/python")

            # Run the MATLAB test.  Note, this doesn't work on Windows
            # since we can't set an environment variable to change the
            # path.
            elseif(${interface} STREQUAL "matlab"
                AND ENABLE_MATLAB_UNIT
                AND NOT WIN32
            )
                # Run the test
                add_test(
                    "Execution_of_matlab_${name}"
                    ${MATLAB_EXECUTABLE}
                    "-nosplash -nodesktop -r"
                    "${name},exit")

                # Make sure Optizelle is available 
                set_tests_properties(
                    "Execution_of_matlab_${name}"
                    PROPERTIES ENVIRONMENT
                    "MATLABPATH=${CMAKE_BINARY_DIR}/src/matlab:${CMAKE_CURRENT_SOURCE_DIR}")
            elseif(${interface} STREQUAL "octave" AND ENABLE_OCTAVE_UNIT)
                # Run the test
                add_test(
                    "Execution_of_octave_${name}"
                    ${OCTAVE_EXECUTABLE}
                    "--eval"
                    "${name},exit")

                # Make sure Optizelle is available 
                set_tests_properties(
                    "Execution_of_octave_${name}"
                    PROPERTIES ENVIRONMENT
                    "OCTAVE_PATH=${CMAKE_BINARY_DIR}/src/octave:${CMAKE_CURRENT_SOURCE_DIR}")
            endif()
        endforeach()
    endif()
endmacro()
