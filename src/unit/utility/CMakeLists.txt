project(utility)

# If any unit test is enabled, we compile this
if(ENABLE_CPP_UNIT)

    # Common include files
    include_directories(${OPTIZELLE_INCLUDE_DIRS})
    include_directories(${JSONCPP_INCLUDE_DIRS})

    # Compile and link the diff utility
    add_executable(diff_restart diff_restart.cpp)
    target_link_libraries(diff_restart
        optizelle_static
        ${JSONCPP_LIBRARIES}
        ${LAPACK_LIBRARIES}
        ${BLAS_LIBRARIES})

endif()
