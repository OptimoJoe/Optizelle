project(simple_constrained)
add_executable(simple_constrained "driver.cpp")
target_link_libraries(simple_constrained
    optizelle_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES driver.cpp simple_constrained.param
        DESTINATION share/optizelle/cpp/examples/simple_constrained)
    install(TARGETS simple_constrained
        DESTINATION share/optizelle/cpp/examples/simple_constrained)
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
    add_test(${param} simple_constrained ${param})
    add_test(${rst}
        "${CMAKE_CURRENT_BINARY_DIR}/../../unit/utility/diff_restart"
        ${rst} simple_constrained.rst)
endforeach()
