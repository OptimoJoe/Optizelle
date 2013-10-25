# Sets up and runs tests on standardard Optizelle unit tests 
macro(add_optizelle_test files executable)
    file(GLOB_RECURSE units ${CMAKE_CURRENT_SOURCE_DIR} ${files})
    list(SORT units)
    list(LENGTH units len1)
    math(EXPR len2 "${len1} - 1")
    foreach(index RANGE ${len2})
        list(GET units ${index} unit)
        add_test("Execution_of_${unit}" ${executable} ${unit})
        add_test("Solution_to_${unit}"
            "${CMAKE_BINARY_DIR}/src/cpp/unit/utility/diff_restart"
            ${unit} solution.json)
    endforeach()
endmacro(add_optizelle_test)
