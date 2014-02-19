# Compiles an Optizelle mex file 
macro(add_optizelle_mex name)
    add_library(${name} SHARED 
        $<TARGET_OBJECTS:optizelle_matlab>
        ${name}.cpp)
    target_link_libraries(${name}
        optizelle_shared
        ${MATLAB_LIBRARIES})
    set_target_properties(${name} PROPERTIES OUTPUT_NAME ${name}_)
    set_target_properties(${name} PROPERTIES PREFIX "") 
    set_target_properties(${name} PROPERTIES SUFFIX .${MATLAB_MEX_EXTENSION}) 
    install(TARGETS ${name} DESTINATION share/optizelle/matlab/optizelle)
    install(FILES ${name}.m DESTINATION share/optizelle/matlab/optizelle)
endmacro()
