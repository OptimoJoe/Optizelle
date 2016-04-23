project(parest)

# Set what interfaces we're using
set(interfaces "matlab" "octave")

# Compile and install the example
set(supporting
    genUnconstrained.m
    genEqualityConstrained.m)
add_example(${PROJECT_NAME} "${interfaces}" "${supporting}")

# Run the unit tests
file(GLOB_RECURSE units ${CMAKE_CURRENT_SOURCE_DIR} "*.json")
add_unit(${PROJECT_NAME} "${interfaces}" "${units}" FALSE)
