project(rosenbrock_custom_VS)
add_executable(rosenbrock_custom_VS "driver.cpp")
target_link_libraries(rosenbrock_custom_VS
    optizelle_static
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})
install(FILES driver.cpp
    DESTINATION share/optizelle/cpp/examples/rosenbrock_custom_VS)
install(TARGETS rosenbrock_custom_VS 
    DESTINATION share/optizelle/cpp/examples/rosenbrock_custom_VS)
