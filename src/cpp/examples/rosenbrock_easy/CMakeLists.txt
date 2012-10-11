project(rosenbrock_easy)
set(rosenbrock_easy_srcs "driver.cpp")
include_directories(${JSONCPP_INCLUDE_DIRS})
add_executable(rosenbrock_easy ${rosenbrock_easy_srcs})
target_link_libraries(rosenbrock_easy
    peopt
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})
install(FILES driver.cpp rosenbrock_easy.peopt
    DESTINATION share/peopt/cpp/examples/rosenbrock_easy)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/rosenbrock_easy
    DESTINATION share/peopt/cpp/examples/rosenbrock_easy)
