project(simple_constrained)
add_executable(simple_constrained "driver.cpp")
target_link_libraries(simple_constrained
    optizelle_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES driver.cpp simple_constrained.json
        DESTINATION share/optizelle/cpp/examples/simple_constrained)
    install(TARGETS simple_constrained
        DESTINATION share/optizelle/cpp/examples/simple_constrained)
endif()

# Run some unit tests
file(GLOB_RECURSE units ${CMAKE_CURRENT_SOURCE_DIR} "unit/*.json")
list(SORT units)
list(LENGTH units len1)
math(EXPR len2 "${len1} - 1")
foreach(index RANGE ${len2})
    list(GET units ${index} unit)
    add_test("Execution_of_${unit}" simple_constrained ${unit})
    add_test("Solution_to_${unit}"
        "${CMAKE_CURRENT_BINARY_DIR}/../../unit/utility/diff_restart"
        ${unit} solution.json)
endforeach()
