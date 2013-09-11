project(simple_sdp_cone)
add_executable(simple_sdp_cone "driver.cpp")
target_link_libraries(simple_sdp_cone
    optizelle_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES driver.cpp simple_sdp_cone.param
        DESTINATION share/optizelle/cpp/examples/simple_sdp_cone)
    install(TARGETS simple_sdp_cone
        DESTINATION share/optizelle/cpp/examples/simple_sdp_cone)
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
    add_test(${param} simple_sdp_cone ${param})
    add_test(${rst}
        "${CMAKE_CURRENT_BINARY_DIR}/../../unit/utility/diff_restart"
        ${rst} simple_sdp_cone.rst)
endforeach()
