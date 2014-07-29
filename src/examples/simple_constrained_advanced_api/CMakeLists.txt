project(simple_constrained_advanced_api)

# Compile and install the example
add_optizelle_example_cpp(${PROJECT_NAME})
add_optizelle_example_python(${PROJECT_NAME})
add_optizelle_example_matlab(${PROJECT_NAME})

# Installs the supporting files 
add_optizelle_example_supporting(${PROJECT_NAME}
    restart.json
    tr_newton.json
    sr1.json)

add_optizelle_example_supporting_restart(${PROJECT_NAME}
    restart/dx.0008.txt
    restart/dx_n.0008.txt
    restart/dx_ncp.0008.txt
    restart/dx_old.0008.txt
    restart/dx_t.0008.txt
    restart/dx_tcp_uncorrected.0008.txt
    restart/dx_t_uncorrected.0008.txt
    restart/dy.0008.txt
    restart/dz.0008.txt
    restart/gpxdxn_p_gx.0008.txt
    restart/gpxdxt.0008.txt
    restart/grad.0008.txt
    restart/grad_old.0008.txt
    restart/g_x.0008.txt
    restart/H_dxn.0008.txt
    restart/H_dxtuncorrected.0008.txt
    restart/h_x.0008.txt
    restart/W_gradpHdxn.0008.txt
    restart/x.0008.txt
    restart/x_old.0008.txt
    restart/y.0008.txt
    restart/z.0008.txt)

# Add some unit tests
add_optizelle_json_test_cpp("*.json" ${PROJECT_NAME}
     "${CMAKE_CURRENT_SOURCE_DIR}/restart.json")
add_optizelle_json_test_python("*.json" ${PROJECT_NAME}
     "${CMAKE_CURRENT_SOURCE_DIR}/restart.json")
add_optizelle_json_test_matlab("*.json" ${PROJECT_NAME}
     "${CMAKE_CURRENT_SOURCE_DIR}/restart.json")
