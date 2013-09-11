project(simple_infeasible_inequality)
add_executable(simple_infeasible_inequality "driver.cpp")
target_link_libraries(simple_infeasible_inequality
    optizelle_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES driver.cpp simple_infeasible_inequality.param
        DESTINATION share/optizelle/cpp/examples/simple_infeasible_inequality)
    install(TARGETS simple_infeasible_inequality
        DESTINATION share/optizelle/cpp/examples/simple_infeasible_inequality)
endif()

# Run some unit tests
file(GLOB_RECURSE params ${CMAKE_CURRENT_SOURCE_DIR} "unit/*.param")
list(SORT params)
file(GLOB_RECURSE rsts ${CMAKE_CURRENT_SOURCE_DIR} "*unit/*.rst")
list(SORT rsts)
list(LENGTH params len1)
math(EXPR len2 "${len1} - 1")
foreach(index RANGE ${len2})
    list(GET params ${index} param)
    list(GET rsts ${index} rst)
    add_test(${param} simple_infeasible_inequality ${param})
    add_test(${rst}
        "${CMAKE_CURRENT_BINARY_DIR}/../../unit/utility/diff_restart"
        ${rst} simple_infeasible_inequality.rst)
endforeach()
