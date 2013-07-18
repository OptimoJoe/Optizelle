project(capture_and_release)

include_directories(${JSONCPP_INCLUDE_DIRS})
include_directories(../utility)

add_executable(unconstrained unconstrained.cpp)
add_executable(equality_constrained equality_constrained.cpp)
add_executable(inequality_constrained inequality_constrained.cpp)
add_executable(constrained constrained.cpp)

target_link_libraries(unconstrained
    peopt_static
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})
target_link_libraries(equality_constrained
    peopt_static
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})
target_link_libraries(inequality_constrained
    peopt_static
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})
target_link_libraries(constrained
    peopt_static
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

add_test(capture_and_release/unconstrained unconstrained)
add_test(capture_and_release/equality_constrained equality_constrained)
add_test(capture_and_release/inequality_constrained inequality_constrained)
add_test(capture_and_release/constrained constrained)
