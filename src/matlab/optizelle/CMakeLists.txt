project(optizelle_matlab)

# Make sure to grab the Matlab includes
include_directories(${MATLAB_INCLUDE_DIRS})

# Compile our core Matlab functionality
set(optizelle_matlab_srcs "optizelle.cpp")
add_library(optizelle_matlab OBJECT ${optizelle_matlab_srcs})

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

# Copy in the rest of the Matlab files for the unit tests. 
add_custom_command(
    TARGET EqualityConstrainedRestartReadRestart POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/src/matlab/optizelle
        $<TARGET_FILE_DIR:EqualityConstrainedRestartReadRestart>)
add_custom_command(
    TARGET EqualityConstrainedRestartReadRestart POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_SOURCE_DIR}/src/matlab/setupOptizelle.m
        $<TARGET_FILE_DIR:EqualityConstrainedRestartReadRestart>/..)

install(FILES
    getStateManipulator.m
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
    DESTINATION share/optizelle/matlab/optizelle)
