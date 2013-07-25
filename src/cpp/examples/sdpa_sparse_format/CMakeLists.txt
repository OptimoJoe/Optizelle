project(sdpa_sparse_format)
add_executable(sdpa_sparse_format "driver.cpp")
target_link_libraries(sdpa_sparse_format
    peopt_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES
        driver.cpp
        sdpa_sparse_format.peopt
        sdpa_sparse_format_phase1.peopt
        lp.dat-s
        DESTINATION share/peopt/cpp/examples/sdpa_sparse_format)
    install(TARGETS sdpa_sparse_format
        DESTINATION share/peopt/cpp/examples/sdpa_sparse_format)
endif()

# Run some unit tests
file(GLOB_RECURSE dats ${CMAKE_CURRENT_SOURCE_DIR} "unit/*.dat-s")
list(SORT dats)
file(GLOB_RECURSE pebases ${CMAKE_CURRENT_SOURCE_DIR} "*unit/*.pebase")
list(SORT pebases)
list(LENGTH dats len1)
math(EXPR len2 "${len1} - 1")
foreach(index RANGE ${len2})
    list(GET dats ${index} dat)
    list(GET pebases ${index} pebase)
    add_test(${dat} sdpa_sparse_format ${dat}
        ${CMAKE_CURRENT_SOURCE_DIR}/sdpa_sparse_format_phase1.peopt
        ${CMAKE_CURRENT_SOURCE_DIR}/sdpa_sparse_format.peopt)
    add_test(${pebase}
        "${CMAKE_CURRENT_BINARY_DIR}/../../unit/utility/diff_restart"
        ${pebase} sdpa_sparse_format.perst)
endforeach()
