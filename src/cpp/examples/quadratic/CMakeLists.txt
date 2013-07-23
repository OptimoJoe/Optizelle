project(quadratic)
include_directories(${JSONCPP_INCLUDE_DIRS})
add_executable(quadratic "driver.cpp")
target_link_libraries(quadratic
    peopt_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES driver.cpp quadratic.peopt
        DESTINATION share/peopt/cpp/examples/quadratic)
    install(TARGETS quadratic
        DESTINATION share/peopt/cpp/examples/quadratic)
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
    add_test(${peopt} quadratic ${peopt})
    add_test(${pebase}
        "${CMAKE_CURRENT_BINARY_DIR}/../../unit/utility/diff_restart"
        ${pebase} quadratic.perst)
endforeach()
