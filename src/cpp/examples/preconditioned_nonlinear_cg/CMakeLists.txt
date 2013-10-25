project(preconditioned_nonlinear_cg)
add_executable(preconditioned_nonlinear_cg "driver.cpp")
target_link_libraries(preconditioned_nonlinear_cg
    optizelle_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES driver.cpp preconditioned_nonlinear_cg.json
        DESTINATION share/optizelle/cpp/examples/preconditioned_nonlinear_cg)
    install(TARGETS preconditioned_nonlinear_cg
        DESTINATION share/optizelle/cpp/examples/preconditioned_nonlinear_cg)
endif()

# Run some unit tests
add_optizelle_test("unit/*.json" preconditioned_nonlinear_cg)
