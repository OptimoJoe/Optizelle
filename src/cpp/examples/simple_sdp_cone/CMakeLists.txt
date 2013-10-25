project(simple_sdp_cone)
add_executable(simple_sdp_cone "driver.cpp")
target_link_libraries(simple_sdp_cone
    optizelle_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES driver.cpp simple_sdp_cone.json
        DESTINATION share/optizelle/cpp/examples/simple_sdp_cone)
    install(TARGETS simple_sdp_cone
        DESTINATION share/optizelle/cpp/examples/simple_sdp_cone)
endif()

# Run some unit tests
file(GLOB_RECURSE units ${CMAKE_CURRENT_SOURCE_DIR} "unit/*.json")
list(SORT units)
list(LENGTH units len1)
math(EXPR len2 "${len1} - 1")
foreach(index RANGE ${len2})
    list(GET units ${index} unit)
    add_test("Exeucution_of_${unit}" simple_sdp_cone ${unit})
    add_test("Solution_to_${unit}"
        "${CMAKE_CURRENT_BINARY_DIR}/../../unit/utility/diff_restart"
        ${unit} solution.json)
endforeach()
