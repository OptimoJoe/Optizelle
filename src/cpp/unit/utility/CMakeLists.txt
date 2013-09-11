project(utility)

add_executable(diff_restart diff_restart.cpp)

target_link_libraries(diff_restart
    optizelle_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})
