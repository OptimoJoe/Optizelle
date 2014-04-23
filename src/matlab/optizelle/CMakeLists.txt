project(optizelle_matlab)

# Make sure to grab the Matlab includes
include_directories(${MATLAB_INCLUDE_DIRS})

# Compile our core Matlab functionality
set(optizelle_matlab_srcs "optizelle.cpp")
add_library(optizelle_matlab SHARED ${optizelle_matlab_srcs})
set_target_properties(optizelle_matlab PROPERTIES OUTPUT_NAME optizelle_matlab)
target_link_libraries(optizelle_matlab 
    optizelle_shared
    ${MATLAB_LIBRARIES})
install(TARGETS optizelle_matlab DESTINATION lib)    

# Compile individual mex files
add_optizelle_mex(UnconstrainedStateCreate)
add_optizelle_mex(UnconstrainedStateReadJson)
add_optizelle_mex(UnconstrainedAlgorithmsGetMin)
add_optizelle_mex(UnconstrainedRestartRelease)
add_optizelle_mex(UnconstrainedRestartCapture)
add_optizelle_mex(UnconstrainedRestartWriteRestart)
add_optizelle_mex(UnconstrainedRestartReadRestart)
add_optizelle_mex(EqualityConstrainedStateCreate)
add_optizelle_mex(EqualityConstrainedStateReadJson)
add_optizelle_mex(EqualityConstrainedAlgorithmsGetMin)
add_optizelle_mex(EqualityConstrainedRestartRelease)
add_optizelle_mex(EqualityConstrainedRestartCapture)
add_optizelle_mex(EqualityConstrainedRestartWriteRestart)
add_optizelle_mex(EqualityConstrainedRestartReadRestart)
add_optizelle_mex(InequalityConstrainedStateCreate)
add_optizelle_mex(InequalityConstrainedStateReadJson)
add_optizelle_mex(InequalityConstrainedAlgorithmsGetMin)
add_optizelle_mex(InequalityConstrainedRestartRelease)
add_optizelle_mex(InequalityConstrainedRestartCapture)
add_optizelle_mex(InequalityConstrainedRestartWriteRestart)
add_optizelle_mex(InequalityConstrainedRestartReadRestart)
add_optizelle_mex(ConstrainedStateCreate)
add_optizelle_mex(ConstrainedStateReadJson)
add_optizelle_mex(ConstrainedAlgorithmsGetMin)
add_optizelle_mex(ConstrainedRestartRelease)
add_optizelle_mex(ConstrainedRestartCapture)
add_optizelle_mex(ConstrainedRestartWriteRestart)
add_optizelle_mex(ConstrainedRestartReadRestart)


# Copy in the rest of the Matlab files for the unit tests. 
if(ENABLE_MATLAB_UNIT)
    add_custom_target(
        matlab_files
        ALL
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${CMAKE_SOURCE_DIR}/src/matlab/optizelle
            ${CMAKE_BINARY_DIR}/src/matlab/optizelle)
    add_custom_target(
        matlab_setup
        ALL
        COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_SOURCE_DIR}/src/matlab/setupOptizelle.m
            ${CMAKE_BINARY_DIR}/src/matlab)
endif()

install(FILES
    serialize.m
    deserialize.m
    checkEnum.m
    checkEnumRange.m
    checkEuclidean.m
    checkFloat.m
    checkFunction.m
    checkMessaging.m
    checkMethod.m
    checkMethods.m
    checkFields.m
    checkNatural.m
    checkOperator.m
    checkScalarValuedFunction.m
    checkStateManipulator.m
    checkString.m
    checkVectorList.m
    checkVectorSpace.m
    checkVectorValuedFunction.m
    checkVectors.m
    checkReals.m
    checkNaturals.m
    checkParams.m
    checkUnconstrainedStateT.m
    checkUnconstrainedFunctionsT.m
    checkEqualityConstrainedStateT.m
    checkEqualityConstrainedFunctionsT.m
    checkInequalityConstrainedStateT.m
    checkInequalityConstrainedFunctionsT.m
    checkConstrainedStateT.m
    checkConstrainedFunctionsT.m
    DESTINATION share/optizelle/matlab/optizelle)
