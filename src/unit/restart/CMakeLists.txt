project(restart)

# Set what interfaces we're using
set(interfaces "cpp" "python" "matlab" "octave")

# Copy in files for our unit tests
is_unit("${interfaces}" enable_unit)
if(enable_unit)
    add_custom_target(
        blank
        ALL
        COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_SOURCE_DIR}/src/unit/restart/blank.json
            ${CMAKE_BINARY_DIR}/src/unit/restart/)
endif()

# Compile and add our unit tests
compile_add_unit(unconstrained "${interfaces}")
compile_add_unit(equality_constrained "${interfaces}")
compile_add_unit(inequality_constrained "${interfaces}")
compile_add_unit(constrained "${interfaces}")
