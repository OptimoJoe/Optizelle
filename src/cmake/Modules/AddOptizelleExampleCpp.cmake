# Compiles an Optizelle C++ example 
macro(add_optizelle_example_cpp name)
    
    # Make sure that examples are enabled 
    if(ENABLE_CPP_EXAMPLES OR ENABLE_CPP_UNIT)

        # Compile and link the example
        add_executable(${name} "${name}.cpp")
        target_link_libraries(${name}
            optizelle_static
            ${JSONCPP_LIBRARIES}
            ${LAPACK_LIBRARIES}
            ${BLAS_LIBRARIES})

        if(ENABLE_CPP_EXAMPLES)
            # Install it to the standard location
            install(FILES "${name}.cpp"
                DESTINATION share/optizelle/cpp/examples/${name})
            install(TARGETS ${name}
                DESTINATION share/optizelle/cpp/examples/${name})
        endif()
    endif()

endmacro()
