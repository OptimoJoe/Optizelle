project(simple_quadratic_cone)
add_executable(simple_quadratic_cone "driver.cpp")
target_link_libraries(simple_quadratic_cone
    optizelle_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES driver.cpp simple_quadratic_cone.json
        DESTINATION share/optizelle/cpp/examples/simple_quadratic_cone)
    install(TARGETS simple_quadratic_cone
        DESTINATION share/optizelle/cpp/examples/simple_quadratic_cone)
endif()

# Run some unit tests
add_optizelle_test("unit/*.json" simple_quadratic_cone)
