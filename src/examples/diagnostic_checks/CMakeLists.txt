project(diagnostic_checks)

# Set what interfaces we're using
set(interfaces "cpp" "python" "matlab" "octave")

# Compile and install the example
add_example(${PROJECT_NAME} "${interfaces}" "")

# Run the unit tests
add_unit(${PROJECT_NAME} "${interfaces}" "" FALSE)
