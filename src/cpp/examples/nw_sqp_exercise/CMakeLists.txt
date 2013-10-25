project(nw_sqp_exercise)
add_executable(nw_sqp_exercise "driver.cpp")
target_link_libraries(nw_sqp_exercise
    optizelle_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES driver.cpp nw_sqp_exercise.param
        DESTINATION share/optizelle/cpp/examples/nw_sqp_exercise)
    install(TARGETS nw_sqp_exercise
        DESTINATION share/optizelle/cpp/examples/nw_sqp_exercise)
endif()

# Run some unit tests
add_optizelle_test("unit/*.json" nw_sqp_exercise)
