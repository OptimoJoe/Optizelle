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
file(GLOB_RECURSE units ${CMAKE_CURRENT_SOURCE_DIR} "unit/*.json")
list(SORT units)
list(LENGTH units len1)
math(EXPR len2 "${len1} - 1")
foreach(index RANGE ${len2})
    list(GET units ${index} unit)
    add_test("Execution_of_${unit}" nw_sqp_exercise ${unit})
    add_test("Solution_to_${unit}"
        "${CMAKE_CURRENT_BINARY_DIR}/../../unit/utility/diff_restart"
        ${unit} solution.json)
endforeach()
