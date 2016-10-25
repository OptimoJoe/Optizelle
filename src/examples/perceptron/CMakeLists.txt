project(perceptron)

# Set what interfaces we're using
set(interfaces "matlab" "octave")

# Compile and install the example
set(supporting
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
add_example(${PROJECT_NAME} "${interfaces}" "${supporting}")

# Run the unit tests
file(COPY ${supporting} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_unit("example_1d" "${interfaces}" "" FALSE)
add_unit("example_1d_extruded" "${interfaces}" "" FALSE)
add_unit("example_2d" "${interfaces}" "" FALSE)
