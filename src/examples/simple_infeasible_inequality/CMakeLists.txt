project(simple_infeasible_inequality)

# Set what interfaces we're using
set(interfaces "cpp" "python" "matlab" "octave")

# Compile and install the example
add_example(${PROJECT_NAME} "${interfaces}" "")

# Run the unit tests
file(GLOB_RECURSE units ${CMAKE_CURRENT_SOURCE_DIR} "*.json")
add_unit(${PROJECT_NAME} "${interfaces}" "${units}" TRUE)
