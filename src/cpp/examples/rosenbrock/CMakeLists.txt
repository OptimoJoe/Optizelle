project(rosenbrock)
set(rosenbrock_srcs "driver.cpp")
add_executable(rosenbrock ${rosenbrock_srcs})
install(FILES driver.cpp
    DESTINATION share/peopt/cpp/examples/rosenbrock)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/rosenbrock 
    DESTINATION share/peopt/cpp/examples/rosenbrock)
