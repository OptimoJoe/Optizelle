project(simple_inequality)
include_directories(${JSONCPP_INCLUDE_DIRS})
add_executable(simple_inequality "driver.cpp")
target_link_libraries(simple_inequality
    peopt_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES driver.cpp simple_inequality.peopt
        DESTINATION share/peopt/cpp/examples/simple_inequality)
    install(TARGETS simple_inequality
        DESTINATION share/peopt/cpp/examples/simple_inequality)
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
    add_test(${peopt} simple_inequality ${peopt})
    add_test(${pebase}
        "${CMAKE_CURRENT_BINARY_DIR}/../../unit/utility/diff_restart"
        ${pebase} simple_inequality.perst)
endforeach()
