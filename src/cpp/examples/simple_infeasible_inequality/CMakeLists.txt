project(simple_infeasible_inequality)
set(simple_infeasible_inequality_srcs "driver.cpp")
include_directories(${JSONCPP_INCLUDE_DIRS})
add_executable(simple_infeasible_inequality
    ${simple_infeasible_inequality_srcs})
target_link_libraries(simple_infeasible_inequality
    peopt_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})
install(FILES driver.cpp simple_infeasible_inequality.peopt
    DESTINATION share/peopt/cpp/examples/simple_infeasible_inequality)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/simple_infeasible_inequality
    DESTINATION share/peopt/cpp/examples/simple_infeasible_inequality)
