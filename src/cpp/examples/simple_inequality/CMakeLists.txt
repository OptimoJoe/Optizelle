project(simple_inequality)
set(simple_inequality_srcs "driver.cpp")
include_directories(${JSONCPP_INCLUDE_DIRS})
add_executable(simple_inequality ${simple_inequality_srcs})
target_link_libraries(simple_inequality
    peopt_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})
install(FILES driver.cpp simple_inequality.peopt
    DESTINATION share/peopt/cpp/examples/simple_inequality)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/simple_inequality
    DESTINATION share/peopt/cpp/examples/simple_inequality)
