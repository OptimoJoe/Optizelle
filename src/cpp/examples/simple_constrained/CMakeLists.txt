project(simple_constrained)
set(simple_constrained_srcs "driver.cpp")
include_directories(${JSONCPP_INCLUDE_DIRS})
add_executable(simple_constrained ${simple_constrained_srcs})
target_link_libraries(simple_constrained
    peopt_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})
install(FILES driver.cpp simple_constrained.peopt
    DESTINATION share/peopt/cpp/examples/simple_constrained)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/simple_constrained
    DESTINATION share/peopt/cpp/examples/simple_constrained)
