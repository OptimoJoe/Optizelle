project(simple_sdp_cone)
include_directories(${JSONCPP_INCLUDE_DIRS})
add_executable(simple_sdp_cone "driver.cpp")
target_link_libraries(simple_sdp_cone
    peopt_static
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})

if(ENABLE_CPP_EXAMPLES)
    install(FILES driver.cpp simple_sdp_cone.peopt
        DESTINATION share/peopt/cpp/examples/simple_sdp_cone)
    install(TARGETS simple_sdp_cone
        DESTINATION share/peopt/cpp/examples/simple_sdp_cone)
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
    add_test(${peopt} simple_sdp_cone ${peopt})
    add_test(${pebase}
        "${CMAKE_CURRENT_BINARY_DIR}/../../unit/utility/diff_restart"
        ${pebase} simple_sdp_cone.perst)
endforeach()
