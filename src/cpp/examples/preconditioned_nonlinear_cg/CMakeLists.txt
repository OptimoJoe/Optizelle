project(preconditioned_nonlinear_cg)
set(preconditioned_nonlinear_cg_srcs "driver.cpp")
include_directories(${JSONCPP_INCLUDE_DIRS})
add_executable(preconditioned_nonlinear_cg ${preconditioned_nonlinear_cg_srcs})
target_link_libraries(preconditioned_nonlinear_cg
    peopt_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})
install(FILES driver.cpp preconditioned_nonlinear_cg.peopt
    DESTINATION share/peopt/cpp/examples/preconditioned_nonlinear_cg)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/preconditioned_nonlinear_cg
    DESTINATION share/peopt/cpp/examples/preconditioned_nonlinear_cg)
