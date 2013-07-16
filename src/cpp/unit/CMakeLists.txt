project(unit)

set(unit_srcs
    "capture_and_release.cpp"
    "linear_cone.cpp"
    "rosenbrock.cpp"
    "quadratic.cpp"
    "sdp_cone.cpp"
    "derivative_check.cpp"
    "quadratic_cone.cpp"
    "linear_algebra.cpp"
    "unit.cpp")
include_directories(${BOOST_UNIT_TEST_INCLUDE_DIRS})
include_directories(${JSONCPP_INCLUDE_DIRS})
add_executable(unit ${unit_srcs})
add_executable(diff_restart diff_restart.cpp)
target_link_libraries(unit
    peopt_static
    ${BOOST_UNIT_TEST_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})
target_link_libraries(diff_restart
    peopt_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})
install(PROGRAMS
    ${CMAKE_CURRENT_BINARY_DIR}/unit
    ${CMAKE_CURRENT_BINARY_DIR}/diff_restart
    DESTINATION share/peopt/cpp/unit)
