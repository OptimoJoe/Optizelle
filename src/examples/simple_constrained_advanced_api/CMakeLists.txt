project(simple_constrained_advanced_api)

# Set what interfaces we're using
set(interfaces "cpp" "python" "matlab" "octave")

# Compile and install the example
set(supporting
    restart.json
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
add_example(${PROJECT_NAME} "${interfaces}" "${supporting}")

# Run the unit tests
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/restart
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(GLOB_RECURSE units ${CMAKE_CURRENT_SOURCE_DIR} "*.json")
list(REMOVE_ITEM units "${CMAKE_CURRENT_SOURCE_DIR}/restart.json") 
add_unit(${PROJECT_NAME} "${interfaces}" "${units}" TRUE
    "${CMAKE_CURRENT_SOURCE_DIR}/restart.json")
