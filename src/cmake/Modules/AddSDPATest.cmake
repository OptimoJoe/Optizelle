# Sets up and runs tests on sparse SDPA formatted files.
macro(add_sdpa_test dats params executable)
    file(GLOB_RECURSE dats ${CMAKE_CURRENT_SOURCE_DIR} ${dats})
    list(SORT dats)
    file(GLOB_RECURSE params ${CMAKE_CURRENT_SOURCE_DIR} ${params})
    list(SORT params)
    list(LENGTH dats len1)
    math(EXPR len2 "${len1} - 1")
    foreach(index RANGE ${len2})
        list(GET dats ${index} dat)
        math(EXPR phase1_idx "2 * ${index}")
        math(EXPR phase2_idx "2 * ${index} + 1")
        list(GET params ${phase1_idx} phase1)
        list(GET params ${phase2_idx} phase2)
        add_test("Execution_of_${dat}"
            ${executable} ${dat} ${phase1} ${phase2})
        add_test("Solution_to_${dat}"
            "${CMAKE_BINARY_DIR}/src/cpp/unit/utility/diff_restart"
            ${phase2} solution.json)
    endforeach()
endmacro(add_sdpa_test)
