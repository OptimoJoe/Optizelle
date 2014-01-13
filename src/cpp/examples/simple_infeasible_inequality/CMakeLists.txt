project(simple_infeasible_inequality)
add_executable(simple_infeasible_inequality "simple_infeasible_inequality.cpp")
target_link_libraries(simple_infeasible_inequality
    optizelle_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES simple_infeasible_inequality.cpp
                  simple_infeasible_inequality.json
        DESTINATION share/optizelle/cpp/examples/simple_infeasible_inequality)
    install(TARGETS simple_infeasible_inequality
        DESTINATION share/optizelle/cpp/examples/simple_infeasible_inequality)
endif()

# Run some unit tests
add_optizelle_test("unit/*.json" simple_infeasible_inequality)
