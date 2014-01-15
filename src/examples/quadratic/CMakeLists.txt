project(quadratic)

# Compile and install the example
add_optizelle_example_cpp(${PROJECT_NAME})
add_optizelle_example_python(${PROJECT_NAME})

# Installs the supporting files 
add_optizelle_example_supporting(${PROJECT_NAME}
    bfgs.json
    fletcher_reeves.json
    hestenes_stiefel.json
    newton_cg.json
    newton_cg_backtracking.json
    polak_ribiere.json
    precon_fletcher_reeves.json
    precon_hestenes_stiefel.json
    precon_polak_ribiere.json
    sr1.json
    tr_newton.json)

# Add some unit tests
add_optizelle_test_cpp("*.json" ${PROJECT_NAME})
add_optizelle_test_python("*.json" ${PROJECT_NAME})
