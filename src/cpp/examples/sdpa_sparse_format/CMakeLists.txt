project(sdpa_sparse_format)

# Compile and install the example
add_optizelle_example_cpp(${PROJECT_NAME})

# Installs the supporting files 
add_optizelle_example_supporting(${PROJECT_NAME}
    lp.dat-s
    lp_phase1.json
    lp_phase2.json)

# Add some unit tests
add_optizelle_sdpa_cpp("*.dat-s" "*.json" ${PROJECT_NAME})
