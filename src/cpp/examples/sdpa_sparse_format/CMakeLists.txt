project(sdpa_sparse_format)
include_directories(${JSONCPP_INCLUDE_DIRS})
add_executable(sdpa_sparse_format "driver.cpp")
target_link_libraries(sdpa_sparse_format
    peopt_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})
install(FILES
    driver.cpp
    sdpa_sparse_format.peopt
    sdpa_sparse_format_phase1.peopt
    lp.dat-s
    DESTINATION share/peopt/cpp/examples/sdpa_sparse_format)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/sdpa_sparse_format
    DESTINATION share/peopt/cpp/examples/sdpa_sparse_format)
