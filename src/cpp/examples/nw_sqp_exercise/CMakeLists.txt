project(nw_sqp_exercise)
add_executable(nw_sqp_exercise "nw_sqp_exercise.cpp")
target_link_libraries(nw_sqp_exercise
    optizelle_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES
        nw_sqp_exercise.cpp
        sr1.json
        tr_newton.json
        DESTINATION share/optizelle/cpp/examples/nw_sqp_exercise)
    install(TARGETS nw_sqp_exercise
        DESTINATION share/optizelle/cpp/examples/nw_sqp_exercise)
endif()

# Run some unit tests
if(ENABLE_CPP_UNIT)
    add_optizelle_test("*.json" nw_sqp_exercise)
endif()
