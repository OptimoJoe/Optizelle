project(derivative_check)

add_executable(scalar_valued scalar_valued.cpp)
add_executable(vector_valued vector_valued.cpp)

target_link_libraries(scalar_valued
    optizelle_static
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})
target_link_libraries(vector_valued
    optizelle_static
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

add_test(derivative_check/scalar_valued scalar_valued)
add_test(derivative_check/vector_valued vector_valued)
