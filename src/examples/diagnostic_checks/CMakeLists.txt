project(diagnostic_checks)

# Compile and install the example
add_optizelle_example_cpp(${PROJECT_NAME})
add_optizelle_example_python(${PROJECT_NAME})
add_optizelle_example_matlab(${PROJECT_NAME})

# Add some unit tests
add_optizelle_test_cpp(${PROJECT_NAME})
add_optizelle_test_python(${PROJECT_NAME})
add_optizelle_test_matlab(${PROJECT_NAME})
