project(simple_inequality)
add_executable(simple_inequality "driver.cpp")
target_link_libraries(simple_inequality
    optizelle_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES driver.cpp simple_inequality.param
        DESTINATION share/optizelle/cpp/examples/simple_inequality)
    install(TARGETS simple_inequality
        DESTINATION share/optizelle/cpp/examples/simple_inequality)
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
    add_test(${param} simple_inequality ${param})
    add_test(${rst}
        "${CMAKE_CURRENT_BINARY_DIR}/../../unit/utility/diff_restart"
        ${rst} simple_inequality.rst)
endforeach()
