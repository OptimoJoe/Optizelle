project(rosenbrock_custom_VS)
set(rosenbrock_srcs "driver.cpp")
add_executable(rosenbrock_custom_VS ${rosenbrock_srcs})
install(FILES driver.cpp
    DESTINATION share/peopt/cpp/examples/rosenbrock_custom_VS)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/rosenbrock_custom_VS 
    DESTINATION share/peopt/cpp/examples/rosenbrock_custom_VS)
