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
        sdpa_sparse_format.param
        sdpa_sparse_format_phase1.param
        lp.dat-s
        DESTINATION share/optizelle/cpp/examples/sdpa_sparse_format)
    install(TARGETS sdpa_sparse_format
        DESTINATION share/optizelle/cpp/examples/sdpa_sparse_format)
endif()

# Run some unit tests
file(GLOB_RECURSE dats ${CMAKE_CURRENT_SOURCE_DIR} "unit/*.dat-s")
list(SORT dats)
file(GLOB_RECURSE rsts ${CMAKE_CURRENT_SOURCE_DIR} "*unit/*.rst")
list(SORT rsts)
list(LENGTH dats len1)
math(EXPR len2 "${len1} - 1")
foreach(index RANGE ${len2})
    list(GET dats ${index} dat)
    list(GET rsts ${index} rst)
    add_test(${dat} sdpa_sparse_format ${dat}
        ${CMAKE_CURRENT_SOURCE_DIR}/sdpa_sparse_format_phase1.param
        ${CMAKE_CURRENT_SOURCE_DIR}/sdpa_sparse_format.param)
    add_test(${rst}
        "${CMAKE_CURRENT_BINARY_DIR}/../../unit/utility/diff_restart"
        ${rst} sdpa_sparse_format.rst)
endforeach()
