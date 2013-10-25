project(rosenbrock)
add_executable(rosenbrock "driver.cpp")
target_link_libraries(rosenbrock
    optizelle_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES driver.cpp rosenbrock.json
        DESTINATION share/optizelle/cpp/examples/rosenbrock)
    install(TARGETS rosenbrock 
        DESTINATION share/optizelle/cpp/examples/rosenbrock)
endif()

# Run some unit tests
add_optizelle_test("unit/*.json" rosenbrock)
