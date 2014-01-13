project(quadratic)
add_executable(quadratic "quadratic.cpp")
target_link_libraries(quadratic
    optizelle_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES
        quadratic.cpp
        bfgs.json
        fletcher_reeves.json
        hestenes_stiefel.json
        newton_cg.json
        newton_cg_backtracking.json
        polak_ribiere.json
        precon_fletcher_reeves.json
        precon_hestenes_stiefel.json
        precon_polak_ribiere.json
        sr1.json
        tr_newton.json
        DESTINATION share/optizelle/cpp/examples/quadratic)
    install(TARGETS quadratic
        DESTINATION share/optizelle/cpp/examples/quadratic)
endif()

# Run some unit tests
add_optizelle_test("*.json" quadratic)
