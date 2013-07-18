project(rosenbrock_custom_VS)
add_executable(rosenbrock_custom_VS "driver.cpp")
install(FILES driver.cpp
    DESTINATION share/peopt/cpp/examples/rosenbrock_custom_VS)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/rosenbrock_custom_VS 
    DESTINATION share/peopt/cpp/examples/rosenbrock_custom_VS)
