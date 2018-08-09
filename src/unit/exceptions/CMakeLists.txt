project(exceptions)

# Set what interfaces we're using
#set(interfaces "cpp" "python" "matlab" "octave")
set(interfaces "cpp" "python" "matlab" "octave")

# Copy in files for our unit tests
is_unit("${interfaces}" enable_unit)
if(enable_unit)
    add_custom_target(
        bad_params
        ALL
        COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_SOURCE_DIR}/src/unit/exceptions/bad_params.json
            ${CMAKE_BINARY_DIR}/src/unit/exceptions/)
endif()

# Copy in files for our unit tests
is_unit("${interfaces}" enable_unit)

# Compile and add our unit tests
compile_add_unit(objective "${interfaces}")
compile_add_unit(read_params "${interfaces}")
