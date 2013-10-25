# Sets up and runs tests on sparse SDPA formatted files.
macro(add_sdpa_test dats params executable)
    file(GLOB_RECURSE dats ${CMAKE_CURRENT_SOURCE_DIR} ${dats})
    list(SORT dats)
    file(GLOB_RECURSE units ${CMAKE_CURRENT_SOURCE_DIR} ${params})
    list(SORT units)
    list(LENGTH dats len1)
    math(EXPR len2 "${len1} - 1")
    foreach(index RANGE ${len2})
        list(GET dats ${index} dat)
        list(GET units ${index} unit)
        add_test("Execution_of_${unit}" ${executable} ${dat}
            ${CMAKE_CURRENT_SOURCE_DIR}/sdpa_sparse_format_phase1.json
            ${CMAKE_CURRENT_SOURCE_DIR}/sdpa_sparse_format.json)
        add_test("Solution_to_${unit}"
            "${CMAKE_BINARY_DIR}/src/cpp/unit/utility/diff_restart"
            ${unit} solution.json)
    endforeach()
endmacro(add_sdpa_test)
