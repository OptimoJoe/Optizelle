project(sdpa_sparse_format)
add_executable(sdpa_sparse_format "driver.cpp")
target_link_libraries(sdpa_sparse_format
    optizelle_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES
        driver.cpp
        sdpa_sparse_format.json
        sdpa_sparse_format_phase1.json
        lp.dat-s
        DESTINATION share/optizelle/cpp/examples/sdpa_sparse_format)
    install(TARGETS sdpa_sparse_format
        DESTINATION share/optizelle/cpp/examples/sdpa_sparse_format)
endif()

# Run some unit tests
file(GLOB_RECURSE dats ${CMAKE_CURRENT_SOURCE_DIR} "unit/*.dat-s")
list(SORT dats)
file(GLOB_RECURSE units ${CMAKE_CURRENT_SOURCE_DIR} "unit/*.json")
list(SORT units)
list(LENGTH dats len1)
math(EXPR len2 "${len1} - 1")
foreach(index RANGE ${len2})
    list(GET dats ${index} dat)
    list(GET units ${index} unit)
    add_test("Execution_of_${unit}" sdpa_sparse_format ${dat}
        ${CMAKE_CURRENT_SOURCE_DIR}/sdpa_sparse_format_phase1.json
        ${CMAKE_CURRENT_SOURCE_DIR}/sdpa_sparse_format.json)
    add_test("Solution_to_${unit}"
        "${CMAKE_CURRENT_BINARY_DIR}/../../unit/utility/diff_restart"
        ${unit} solution.json)
endforeach()
