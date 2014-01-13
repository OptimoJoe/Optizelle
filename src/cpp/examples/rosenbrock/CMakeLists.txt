project(rosenbrock)
add_executable(rosenbrock "rosenbrock.cpp")
target_link_libraries(rosenbrock
    optizelle_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES
        rosenbrock.cpp
        bfgs.json
        newton_cg.json
        newton_cg_backtracking.json
        precon_fletcher_reeves.json
        precon_hestenes_stiefel.json
        precon_polak_ribiere.json
        sr1.json
        tr_newton.json
        DESTINATION share/optizelle/cpp/examples/rosenbrock)
    install(TARGETS rosenbrock 
        DESTINATION share/optizelle/cpp/examples/rosenbrock)
endif()

# Run some unit tests
add_optizelle_test("*.json" rosenbrock)
