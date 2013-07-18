project(rosenbrock)
include_directories(${JSONCPP_INCLUDE_DIRS})
add_executable(rosenbrock "driver.cpp")
target_link_libraries(rosenbrock
    peopt_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES driver.cpp rosenbrock.peopt
        DESTINATION share/peopt/cpp/examples/rosenbrock)
    install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/rosenbrock
        DESTINATION share/peopt/cpp/examples/rosenbrock)
endif()

# Run some unit tests
file(GLOB_RECURSE peopts ${CMAKE_CURRENT_SOURCE_DIR} "unit/*.peopt")
list(SORT peopts)
file(GLOB_RECURSE pebases ${CMAKE_CURRENT_SOURCE_DIR} "*unit/*.pebase")
list(SORT pebases)
list(LENGTH peopts len1)
math(EXPR len2 "${len1} - 1")
foreach(index RANGE ${len2})
    list(GET peopts ${index} peopt)
    list(GET pebases ${index} pebase)
    add_test(${peopt} rosenbrock ${peopt})
    add_test(${pebase}
        "${CMAKE_CURRENT_BINARY_DIR}/../../unit/utility/diff_restart"
        ${pebase} rosenbrock.perst)
endforeach()
