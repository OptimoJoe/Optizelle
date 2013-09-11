project(diagnostic_checks)
set(diagnostic_checks_srcs "driver.cpp")
add_executable(diagnostic_checks "driver.cpp")
target_link_libraries(diagnostic_checks optizelle_static)
install(FILES driver.cpp
    DESTINATION share/optizelle/cpp/examples/diagnostic_checks)
install(TARGETS diagnostic_checks 
    DESTINATION share/optizelle/cpp/examples/diagnostic_checks)
