project(simple_sdp_cone)
set(simple_sdp_cone_srcs "driver.cpp")
include_directories(${JSONCPP_INCLUDE_DIRS})
add_executable(simple_sdp_cone ${simple_sdp_cone_srcs})
target_link_libraries(simple_sdp_cone
    peopt_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})
install(FILES driver.cpp simple_sdp_cone.peopt
    DESTINATION share/peopt/cpp/examples/simple_sdp_cone)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/simple_sdp_cone
    DESTINATION share/peopt/cpp/examples/simple_sdp_cone)
