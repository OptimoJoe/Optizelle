project(simple_equality)
add_executable(simple_equality "driver.cpp")
target_link_libraries(simple_equality
    optizelle_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES driver.cpp simple_equality.param
        DESTINATION share/optizelle/cpp/examples/simple_equality)
    install(TARGETS simple_equality
        DESTINATION share/optizelle/cpp/examples/simple_equality)
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
    add_test(${param} simple_equality ${param})
    add_test(${rst}
        "${CMAKE_CURRENT_BINARY_DIR}/../../unit/utility/diff_restart"
        ${rst} simple_equality.rst)
endforeach()
