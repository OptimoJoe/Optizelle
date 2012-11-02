project(quadratic)
set(quadratic_srcs "driver.cpp")
include_directories(${JSONCPP_INCLUDE_DIRS})
add_executable(quadratic ${quadratic_srcs})
target_link_libraries(quadratic
    peopt_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})
install(FILES driver.cpp quadratic.peopt
    DESTINATION share/peopt/cpp/examples/quadratic)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/quadratic
    DESTINATION share/peopt/cpp/examples/quadratic)
