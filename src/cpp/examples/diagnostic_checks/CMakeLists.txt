project(diagnostic_checks)
set(diagnostic_checks_srcs "diagnostic_checks.cpp")
add_executable(diagnostic_checks "diagnostic_checks.cpp")
target_link_libraries(diagnostic_checks optizelle_static)
install(FILES diagnostic_checks.cpp
    DESTINATION share/optizelle/cpp/examples/diagnostic_checks)
install(TARGETS diagnostic_checks 
    DESTINATION share/optizelle/cpp/examples/diagnostic_checks)
