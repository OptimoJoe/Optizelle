project(rosenbrock_custom_VS)

# Compile and install the example
add_optizelle_example_cpp(${PROJECT_NAME})
add_optizelle_example_supporting(${PROJECT_NAME}
    rosenbrock_custom_VS.py)
