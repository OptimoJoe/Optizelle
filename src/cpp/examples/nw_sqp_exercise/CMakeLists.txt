project(nw_sqp_exercise)
add_executable(nw_sqp_exercise "driver.cpp")
target_link_libraries(nw_sqp_exercise
    peopt_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES driver.cpp nw_sqp_exercise.peopt
        DESTINATION share/peopt/cpp/examples/nw_sqp_exercise)
    install(TARGETS nw_sqp_exercise
        DESTINATION share/peopt/cpp/examples/nw_sqp_exercise)
endif()

# Run some unit tests
file(GLOB_RECURSE peopts ${CMAKE_CURRENT_SOURCE_DIR} "unit/*.peopt")
list(SORT peopts)
file(GLOB_RECURSE pebases ${CMAKE_CURRENT_SOURCE_DIR} "*unit/*.pebase")
list(SORT pebases)
list(LENGTH peopts len1)
math(EXPR len2 "${len1} - 1")
foreach(index RANGE ${len2})
    list(GET peopts ${index} peopt)
    list(GET pebases ${index} pebase)
    add_test(${peopt} nw_sqp_exercise ${peopt})
    add_test(${pebase}
        "${CMAKE_CURRENT_BINARY_DIR}/../../unit/utility/diff_restart"
        ${pebase} nw_sqp_exercise.perst)
endforeach()
