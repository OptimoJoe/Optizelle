project(simple_equality)
add_executable(simple_equality "simple_equality.cpp")
target_link_libraries(simple_equality
    optizelle_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES simple_equality.cpp simple_equality.json
        DESTINATION share/optizelle/cpp/examples/simple_equality)
    install(TARGETS simple_equality
        DESTINATION share/optizelle/cpp/examples/simple_equality)
endif()

# Run some unit tests
add_optizelle_test("unit/*.json" simple_equality)
