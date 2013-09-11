project(optizelle_python)

set(optizelle_python_srcs "optizelle.cpp")
include_directories(${PYTHON_INCLUDE_DIRS})
add_library(optizelle_python SHARED
    $<TARGET_OBJECTS:optizelle_cpp>
    ${optizelle_python_srcs})
target_link_libraries(optizelle_python
    optizelle_static
    ${PYTHON_LIBRARIES}
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})
set_target_properties(optizelle_python PROPERTIES OUTPUT_NAME optizelle_python)
set_target_properties(optizelle_python PROPERTIES PREFIX "") 
set_target_properties(optizelle_python PROPERTIES SUFFIX .pylib) 
install(FILES Optizelle.py DESTINATION share/optizelle/python)
install(TARGETS optizelle_python DESTINATION share/optizelle/python)
