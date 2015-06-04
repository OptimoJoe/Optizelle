project(perceptron)

# Installs the supporting files 
add_optizelle_example_supporting(${PROJECT_NAME}
    diagnostic_checks_auxiliary.m
    diagnostic_checks_objective.m
    example_1d_extruded.json
    example_1d_extruded.m
    example_1d.json
    example_1d.m
    example_2d.json
    example_2d.m
    generate_hyperbolic.m
    generate_interpolant.m
    generate_lenses.m
    generate_logistic.m
    generate_objective.m
    generate_parametrization.m
    generate_scaling.m
    generate_scaling_operator.m)

# Copy in the rest of the files for the unit tests. 
if(ENABLE_MATLAB_UNIT)
    file(GLOB perceptron_json
        ${CMAKE_SOURCE_DIR}/src/examples/perceptron/*.json)
    file(COPY
        ${perceptron_json}
        DESTINATION ${CMAKE_BINARY_DIR}/src/examples/perceptron)
endif()

# Add some unit tests
add_optizelle_test_matlab("example_1d" , TRUE)
add_optizelle_test_matlab("example_1d_extruded" , TRUE)
add_optizelle_test_matlab("example_2d" , TRUE)
