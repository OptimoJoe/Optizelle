project(computation_caching)

# Set what interfaces we're using
set(interfaces "matlab" "octave")

# Compile and install the example
set(supporting
    driver_full_space.m
    driver_reduced_space.m
    full_space.json
    generate_convection_bd.m
    generate_convection.m
    generate_diffusion_bd.m
    generate_diffusion.m
    generate_full_space.m
    generate_params.m
    generate_reduced_space.m
    reduced_space.json)
add_example(${PROJECT_NAME} "${interfaces}" "${supporting}")

# Run the unit tests
file(COPY ${supporting} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_unit("driver_reduced_space" "${interfaces}" "" FALSE)
add_unit("driver_full_space" "${interfaces}" "" FALSE)
