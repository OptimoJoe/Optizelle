project(preconditioned_nonlinear_cg)
add_executable(preconditioned_nonlinear_cg "driver.cpp")
target_link_libraries(preconditioned_nonlinear_cg
    optizelle_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES driver.cpp preconditioned_nonlinear_cg.param
        DESTINATION share/optizelle/cpp/examples/preconditioned_nonlinear_cg)
    install(TARGETS preconditioned_nonlinear_cg
        DESTINATION share/optizelle/cpp/examples/preconditioned_nonlinear_cg)
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
    add_test(${param} preconditioned_nonlinear_cg ${param})
    add_test(${rst}
        "${CMAKE_CURRENT_BINARY_DIR}/../../unit/utility/diff_restart"
        ${rst} preconditioned_nonlinear_cg.rst)
endforeach()
