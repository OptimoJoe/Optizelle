project(nw_sqp_exercise)
add_executable(nw_sqp_exercise "driver.cpp")
target_link_libraries(nw_sqp_exercise
    optizelle_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES driver.cpp nw_sqp_exercise.param
        DESTINATION share/optizelle/cpp/examples/nw_sqp_exercise)
    install(TARGETS nw_sqp_exercise
        DESTINATION share/optizelle/cpp/examples/nw_sqp_exercise)
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
    add_test(${param} nw_sqp_exercise ${param})
    add_test(${rst}
        "${CMAKE_CURRENT_BINARY_DIR}/../../unit/utility/diff_restart"
        ${rst} nw_sqp_exercise.rst)
endforeach()
