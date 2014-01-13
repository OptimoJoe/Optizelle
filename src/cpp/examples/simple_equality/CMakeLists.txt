project(simple_equality)
add_executable(simple_equality "simple_equality.cpp")
target_link_libraries(simple_equality
    optizelle_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES
        simple_equality.cpp
        sr1.json
        tr_newton.json
        DESTINATION share/optizelle/cpp/examples/simple_equality)
    install(TARGETS simple_equality
        DESTINATION share/optizelle/cpp/examples/simple_equality)
endif()

# Run some unit tests
if(ENABLE_CPP_UNIT)
    add_optizelle_test("*.json" simple_equality)
endif()
