project(sdpa_sparse_format)

# Set what interfaces we're using
set(interfaces "cpp")

# Compile and install the example
file(GLOB_RECURSE supporting ${CMAKE_CURRENT_SOURCE_DIR} "*.dat-s" "*.json")
add_example(${PROJECT_NAME} "${interfaces}" "${supporting}")

# Run the unit tests
if(ENABLE_CPP_UNIT)
    # Grab the problem setups
    file(GLOB_RECURSE units ${CMAKE_CURRENT_SOURCE_DIR} "*.dat-s")

    # Run the tests
    foreach(unit0 ${units})
        string(REGEX REPLACE "\.dat-s" "" unit ${unit0})
        get_filename_component(uname ${unit} NAME_WE)
        add_test("Execution_of_cpp_${PROJECT_NAME}_${uname}"
            ${PROJECT_NAME}
            "${unit}.dat-s"
            "${unit}_phase1.json"
            "${unit}_phase2.json")
    endforeach()
endif()
