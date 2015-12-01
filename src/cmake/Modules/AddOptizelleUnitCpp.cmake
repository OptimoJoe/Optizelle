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
            optizelle_shared)
            
        # Grab the absolute path of the test
        get_filename_component(fullname ${name} ABSOLUTE)

        # Add the test
        add_test("Execution_of_cpp_${fullname}" ${name})
    endif()

endmacro()
