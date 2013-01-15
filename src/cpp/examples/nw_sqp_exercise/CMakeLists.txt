project(nw_sqp_exercise)
set(nw_sqp_exercise_srcs "driver.cpp")
include_directories(${JSONCPP_INCLUDE_DIRS})
add_executable(nw_sqp_exercise ${nw_sqp_exercise_srcs})
target_link_libraries(nw_sqp_exercise
    peopt_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})
install(FILES driver.cpp nw_sqp_exercise.peopt
    DESTINATION share/peopt/cpp/examples/nw_sqp_exercise)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/nw_sqp_exercise
    DESTINATION share/peopt/cpp/examples/nw_sqp_exercise)
