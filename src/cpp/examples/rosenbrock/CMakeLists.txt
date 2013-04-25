project(rosenbrock)
set(rosenbrock_srcs "driver.cpp")
include_directories(${JSONCPP_INCLUDE_DIRS})
add_executable(rosenbrock ${rosenbrock_srcs})
target_link_libraries(rosenbrock
    peopt_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})
install(FILES driver.cpp rosenbrock.peopt
    DESTINATION share/peopt/cpp/examples/rosenbrock)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/rosenbrock
    DESTINATION share/peopt/cpp/examples/rosenbrock)
