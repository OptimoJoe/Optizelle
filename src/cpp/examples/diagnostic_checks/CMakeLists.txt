project(diagnostic_checks)
set(diagnostic_checks_srcs "driver.cpp")
add_executable(diagnostic_checks "driver.cpp")
target_link_libraries(diagnostic_checks peopt_static)
install(FILES driver.cpp
    DESTINATION share/peopt/cpp/examples/diagnostic_checks)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/diagnostic_checks 
    DESTINATION share/peopt/cpp/examples/diagnostic_checks)
