project(parest)

# Compile and install the example
add_optizelle_example_matlab(${PROJECT_NAME})

# Installs the supporting files 
add_optizelle_example_supporting(${PROJECT_NAME}
    genUnconstrained.m
    genEqualityConstrained.m
    tr_newton.json)
