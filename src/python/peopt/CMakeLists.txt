project(peopt_python)

set(peopt_python_srcs "peopt.cpp")
include_directories(${PYTHON_INCLUDE_DIRS})
add_library(peopt_python SHARED
    $<TARGET_OBJECTS:peopt_cpp>
    ${peopt_python_srcs})
target_link_libraries(peopt_python
    peopt_static
    ${PYTHON_LIBRARIES}
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})
set_target_properties(peopt_python PROPERTIES OUTPUT_NAME peopt_python)
set_target_properties(peopt_python PROPERTIES PREFIX "") 
set_target_properties(peopt_python PROPERTIES SUFFIX .pylib) 
install(FILES peopt.py DESTINATION share/peopt/python)
install(TARGETS peopt_python DESTINATION share/peopt/python)
