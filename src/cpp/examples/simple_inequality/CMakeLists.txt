project(simple_inequality)
add_executable(simple_inequality "driver.cpp")
target_link_libraries(simple_inequality
    optizelle_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES driver.cpp simple_inequality.json
        DESTINATION share/optizelle/cpp/examples/simple_inequality)
    install(TARGETS simple_inequality
        DESTINATION share/optizelle/cpp/examples/simple_inequality)
endif()

# Run some unit tests
add_optizelle_test("unit/*.json" simple_inequality)
