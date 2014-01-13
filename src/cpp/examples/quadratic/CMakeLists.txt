project(quadratic)
add_executable(quadratic "quadratic.cpp")
target_link_libraries(quadratic
    optizelle_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES quadratic.cpp quadratic.json
        DESTINATION share/optizelle/cpp/examples/quadratic)
    install(TARGETS quadratic
        DESTINATION share/optizelle/cpp/examples/quadratic)
endif()

# Run some unit tests
add_optizelle_test("unit/*.json" quadratic)
