project(derivative_check)

include_directories(${JSONCPP_INCLUDE_DIRS})
include_directories(../utility)

add_executable(scalar_valued scalar_valued.cpp)
add_executable(vector_valued vector_valued.cpp)

target_link_libraries(scalar_valued
    peopt_static
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})
target_link_libraries(vector_valued
    peopt_static
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

add_test(derivative_check/scalar_valued scalar_valued)
add_test(derivative_check/vector_valued vector_valued)
