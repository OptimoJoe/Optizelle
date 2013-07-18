project(utility)

include_directories(${JSONCPP_INCLUDE_DIRS})

add_executable(diff_restart diff_restart.cpp)

target_link_libraries(diff_restart
    peopt_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})
