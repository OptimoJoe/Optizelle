cmake_minimum_required(VERSION 2.8.9)
project(unit)

set(unit_srcs
    "capture_and_release.cpp"
    "linear_cone.cpp"
    "rosenbrock.cpp"
    "sdp_cone.cpp"
    "derivative_check.cpp"
    "quadratic_cone.cpp"
    "unit.cpp")
include_directories(${BOOST_UNIT_TEST_INCLUDE_DIRS})
add_executable(unit ${unit_srcs})
target_link_libraries(unit
    peopt_static
    ${BOOST_UNIT_TEST_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/unit
    DESTINATION share/peopt/cpp/unit)
