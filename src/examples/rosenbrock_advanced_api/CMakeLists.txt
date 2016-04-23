project(rosenbrock_advanced_api)

# Set what interfaces we're using
set(interfaces "cpp" "python" "matlab" "octave")

# Compile and install the example
add_example(${PROJECT_NAME} "${interfaces}" "restart.json")

# Run the unit tests
file(GLOB_RECURSE units ${CMAKE_CURRENT_SOURCE_DIR} "*.json")
list(REMOVE_ITEM units "${CMAKE_CURRENT_SOURCE_DIR}/restart.json") 
add_unit(${PROJECT_NAME} "${interfaces}" "${units}" TRUE
    "${CMAKE_CURRENT_SOURCE_DIR}/restart.json")
