project(rosenbrock)
add_executable(rosenbrock "driver.cpp")
target_link_libraries(rosenbrock
    optizelle_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES driver.cpp rosenbrock.json
        DESTINATION share/optizelle/cpp/examples/rosenbrock)
    install(TARGETS rosenbrock 
        DESTINATION share/optizelle/cpp/examples/rosenbrock)
endif()

# Run some unit tests
file(GLOB_RECURSE units ${CMAKE_CURRENT_SOURCE_DIR} "unit/*.json")
list(SORT units)
list(LENGTH units len1)
math(EXPR len2 "${len1} - 1")
foreach(index RANGE ${len2})
    list(GET units ${index} unit)
    add_test("Execution_of_${unit}" rosenbrock ${unit})
    add_test("Solution_to_${unit}"
        "${CMAKE_CURRENT_BINARY_DIR}/../../unit/utility/diff_restart"
        ${unit} solution.json)
endforeach()
