project(simple_inequality)
add_executable(simple_inequality "simple_inequality.cpp")
target_link_libraries(simple_inequality
    optizelle_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES
        simple_inequality.cpp
        newton_cg.json
        newton_cg_backtracking.json
        sr1.json
        tr_newton.json
        tr_newton_predictor_corrector.json
        DESTINATION share/optizelle/cpp/examples/simple_inequality)
    install(TARGETS simple_inequality
        DESTINATION share/optizelle/cpp/examples/simple_inequality)
endif()

# Run some unit tests
add_optizelle_test("*.json" simple_inequality)
