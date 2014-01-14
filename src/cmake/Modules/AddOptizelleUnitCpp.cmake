# Compiles an Optizelle C++ unit test 
macro(add_optizelle_unit_cpp name)
    
    # Make sure that unit tests are enabled 
    if(ENABLE_CPP_UNIT)

        # Set common includes
        include_directories(${OPTIZELLE_INCLUDE_DIRS})
        include_directories(${JSONCPP_INCLUDE_DIRS})

        # Compile and link the example
        add_executable(${name} "${name}.cpp")
        target_link_libraries(${name}
            optizelle_static
            ${JSONCPP_LIBRARIES}
            ${LAPACK_LIBRARIES}
            ${BLAS_LIBRARIES})
            
        add_test("Execution_of_${name}" ${name})
    endif()

endmacro()
