project(simple_constrained)

# Compile and install the example
add_optizelle_example_cpp(${PROJECT_NAME})

# Installs the supporting files 
add_optizelle_example_supporting(${PROJECT_NAME}
    tr_newton.json
    sr1.json)

# Add some unit tests
add_optizelle_test_cpp("*.json" ${PROJECT_NAME})
