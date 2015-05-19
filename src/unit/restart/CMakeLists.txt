project(restart)

# Copy in files for our unit tests 
if(ENABLE_CPP_UNIT OR ENABLE_PYTHON_UNIT OR ENABLE_MATLAB_UNIT)
    add_custom_target(
        blank
        ALL
        COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_SOURCE_DIR}/src/unit/restart/blank.json
            ${CMAKE_BINARY_DIR}/src/unit/restart/)
endif()

add_optizelle_unit_cpp(unconstrained)
add_optizelle_test_python(unconstrained)
add_optizelle_test_matlab(unconstrained)

add_optizelle_unit_cpp(equality_constrained)
add_optizelle_test_python(equality_constrained)
add_optizelle_test_matlab(equality_constrained)

add_optizelle_unit_cpp(inequality_constrained)
add_optizelle_test_python(inequality_constrained)
add_optizelle_test_matlab(inequality_constrained)

add_optizelle_unit_cpp(constrained)
add_optizelle_test_python(constrained)
add_optizelle_test_matlab(constrained)
