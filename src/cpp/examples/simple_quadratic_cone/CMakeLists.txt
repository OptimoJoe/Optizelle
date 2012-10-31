project(simple_quadratic_cone)
set(simple_quadratic_cone_srcs "driver.cpp")
include_directories(${JSONCPP_INCLUDE_DIRS})
add_executable(simple_quadratic_cone ${simple_quadratic_cone_srcs})
target_link_libraries(simple_quadratic_cone
    peopt_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})
install(FILES driver.cpp simple_quadratic_cone.peopt
    DESTINATION share/peopt/cpp/examples/simple_quadratic_cone)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/simple_quadratic_cone
    DESTINATION share/peopt/cpp/examples/simple_quadratic_cone)
