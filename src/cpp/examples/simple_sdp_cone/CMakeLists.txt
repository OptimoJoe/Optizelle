project(simple_sdp_cone)
add_executable(simple_sdp_cone "simple_sdp_cone.cpp")
target_link_libraries(simple_sdp_cone
    optizelle_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES simple_sdp_cone.cpp simple_sdp_cone.json
        DESTINATION share/optizelle/cpp/examples/simple_sdp_cone)
    install(TARGETS simple_sdp_cone
        DESTINATION share/optizelle/cpp/examples/simple_sdp_cone)
endif()

# Run some unit tests
add_optizelle_test("unit/*.json" simple_sdp_cone)
