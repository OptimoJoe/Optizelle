project(simple_constrained)
add_executable(simple_constrained "simple_constrained.cpp")
target_link_libraries(simple_constrained
    optizelle_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES
        simple_constrained.cpp
        tr_newton.json
        sr1.json
        DESTINATION share/optizelle/cpp/examples/simple_constrained)
    install(TARGETS simple_constrained
        DESTINATION share/optizelle/cpp/examples/simple_constrained)
endif()

# Run some unit tests
add_optizelle_test("*.json" simple_constrained)
