project(simple_equality)
set(simple_equality_srcs "driver.cpp")
include_directories(${JSONCPP_INCLUDE_DIRS})
add_executable(simple_equality ${simple_equality_srcs})
target_link_libraries(simple_equality
    peopt_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})
install(FILES driver.cpp simple_equality.peopt
    DESTINATION share/peopt/cpp/examples/simple_equality)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/simple_equality
    DESTINATION share/peopt/cpp/examples/simple_equality)
