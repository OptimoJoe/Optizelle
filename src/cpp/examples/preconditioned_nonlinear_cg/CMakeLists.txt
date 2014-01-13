project(preconditioned_nonlinear_cg)
add_executable(preconditioned_nonlinear_cg "preconditioned_nonlinear_cg.cpp")
target_link_libraries(preconditioned_nonlinear_cg
    optizelle_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES
        preconditioned_nonlinear_cg.cpp
        fletcher_reeves.json
        tr_newton.json
        DESTINATION share/optizelle/cpp/examples/preconditioned_nonlinear_cg)
    install(TARGETS preconditioned_nonlinear_cg
        DESTINATION share/optizelle/cpp/examples/preconditioned_nonlinear_cg)
endif()

# Run some unit tests
if(ENABLE_CPP_UNIT)
    add_optizelle_test("*.json" preconditioned_nonlinear_cg)
endif()
