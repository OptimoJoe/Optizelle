# Many of these examples are also unit test
if(ENABLE_CPP_EXAMPLES OR ENABLE_CPP_UNIT)
    add_subdirectory(rosenbrock)
    add_subdirectory(quadratic)
    add_subdirectory(simple_inequality)
    add_subdirectory(simple_infeasible_inequality)
    add_subdirectory(simple_quadratic_cone)
    add_subdirectory(simple_sdp_cone)
    add_subdirectory(simple_equality)
    add_subdirectory(nw_sqp_exercise)
    add_subdirectory(simple_constrained)
    add_subdirectory(sdpa_sparse_format)
endif()

# Though, some examples are just examples
if(ENABLE_CPP_EXAMPLES) 
    add_subdirectory(diagnostic_checks)
    add_subdirectory(rosenbrock_custom_VS)
endif()
