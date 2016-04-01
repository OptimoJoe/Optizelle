project(computation_caching)

# Installs the supporting files 
add_optizelle_example_supporting(${PROJECT_NAME}
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

# Copy in the rest of the files for the unit tests. 
if(ENABLE_MATLAB_UNIT)
    file(GLOB computation_caching_json
        ${CMAKE_SOURCE_DIR}/src/examples/computation_caching/*.json)
    file(COPY
        ${computation_caching_json}
        DESTINATION ${CMAKE_BINARY_DIR}/src/examples/computation_caching)
endif()

# Add some unit tests
add_optizelle_test_matlab("driver_reduced_space" , TRUE)
add_optizelle_test_matlab("driver_full_space" , TRUE)
