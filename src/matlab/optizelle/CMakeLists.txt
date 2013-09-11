project(optizelle_matlab)

set(optizelle_matlab_srcs "optizelle.cpp")
include_directories(${MATLAB_INCLUDE_DIRS})
add_library(optizelle_matlab SHARED ${optizelle_matlab_srcs})
target_link_libraries(optizelle_matlab
    optizelle_static
    ${MATLAB_LIBRARIES}
    ${JSONCPP_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES})
set_target_properties(optizelle_matlab PROPERTIES OUTPUT_NAME Optizelle)
set_target_properties(optizelle_matlab PROPERTIES PREFIX "") 
set_target_properties(optizelle_matlab PROPERTIES SUFFIX .${MATLAB_MEX_EXTENSION}) 
install(TARGETS optizelle_matlab DESTINATION share/optizelle/matlab)
