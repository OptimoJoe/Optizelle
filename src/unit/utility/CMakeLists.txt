project(utility)

# Set what interfaces we're using
set(interfaces "cpp" "python" "matlab" "octave")

# If any unit test is enabled, we compile this
is_unit("${interfaces}" enable_unit)
if(enable_unit)

    # Common include files
    include_directories(${OPTIZELLE_INCLUDE_DIRS})
    include_directories(${JSONCPP_INCLUDE_DIRS})

    # Compile and link the diff utility
    add_executable(diff_restart diff_restart.cpp)
    target_link_libraries(diff_restart optizelle_shared)
endif()
