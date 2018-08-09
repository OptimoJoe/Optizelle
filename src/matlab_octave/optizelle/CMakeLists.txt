# Compiles a mex file
macro(add_mex name type include extension)
    add_library(${name}_${type} SHARED ${name}.cpp)
    target_include_directories(${name}_${type} PRIVATE ${include})
    target_link_libraries(${name}_${type} optizelle_${type})
    set_target_properties(${name}_${type} PROPERTIES OUTPUT_NAME ${name}_)
    set_target_properties(${name}_${type} PROPERTIES PREFIX "")
    set_target_properties(${name}_${type} PROPERTIES SUFFIX .${extension})
    if(APPLE)
        set_target_properties(${name}_${type} PROPERTIES INSTALL_RPATH "@loader_path/../../../../lib;@loader_path/../../thirdparty/lib")
    elseif(UNIX)
        set_target_properties(${name}_${type} PROPERTIES INSTALL_RPATH "\$ORIGIN/../../../../lib;\$ORIGIN/../../thirdparty/lib")
    endif()
    install(TARGETS ${name}_${type}
        DESTINATION share/optizelle/${type}/optizelle)
    install(FILES ${name}.m DESTINATION share/optizelle/${type}/optizelle)
endmacro()

# Compiles the MATLAB/Octave interface
macro(compile_interface type include libs extension unit)

    # Compile our core functionality
    set(optizelle_srcs "optizelle.cpp")
    add_library(optizelle_${type} SHARED ${optizelle_srcs})
    target_include_directories(optizelle_${type} PRIVATE ${include})
    set_target_properties(optizelle_${type}
        PROPERTIES OUTPUT_NAME "optizelle_${type}")
    target_link_libraries(optizelle_${type} optizelle_shared ${libs})
    install(TARGETS optizelle_${type} DESTINATION lib)

    # Compile individual mex files
    add_mex(UnconstrainedStateCreate ${type} ${include} ${extension})
    add_mex(UnconstrainedStateReadJson ${type} ${include} ${extension})
    add_mex(UnconstrainedAlgorithmsGetMin ${type} ${include} ${extension})
    add_mex(UnconstrainedRestartRelease ${type} ${include} ${extension})
    add_mex(UnconstrainedRestartCapture ${type} ${include} ${extension})
    add_mex(UnconstrainedRestartWriteRestart ${type} ${include} ${extension})
    add_mex(UnconstrainedRestartReadRestart ${type} ${include} ${extension})
    add_mex(EqualityConstrainedStateCreate ${type} ${include} ${extension})
    add_mex(EqualityConstrainedStateReadJson ${type} ${include} ${extension})
    add_mex(EqualityConstrainedAlgorithmsGetMin ${type} ${include} ${extension})
    add_mex(EqualityConstrainedRestartRelease ${type} ${include} ${extension})
    add_mex(EqualityConstrainedRestartCapture ${type} ${include} ${extension})
    add_mex(EqualityConstrainedRestartWriteRestart ${type} ${include}
        ${extension})
    add_mex(EqualityConstrainedRestartReadRestart ${type} ${include}
        ${extension})
    add_mex(InequalityConstrainedStateCreate ${type} ${include} ${extension})
    add_mex(InequalityConstrainedStateReadJson ${type} ${include} ${extension})
    add_mex(InequalityConstrainedAlgorithmsGetMin ${type} ${include}
        ${extension})
    add_mex(InequalityConstrainedRestartRelease ${type} ${include} ${extension})
    add_mex(InequalityConstrainedRestartCapture ${type} ${include} ${extension})
    add_mex(InequalityConstrainedRestartWriteRestart ${type} ${include}
        ${extension})
    add_mex(InequalityConstrainedRestartReadRestart ${type} ${include}
        ${extension})
    add_mex(ConstrainedStateCreate ${type} ${include} ${extension})
    add_mex(ConstrainedStateReadJson ${type} ${include} ${extension})
    add_mex(ConstrainedAlgorithmsGetMin ${type} ${include} ${extension})
    add_mex(ConstrainedRestartRelease ${type} ${include} ${extension})
    add_mex(ConstrainedRestartCapture ${type} ${include} ${extension})
    add_mex(ConstrainedRestartWriteRestart ${type} ${include} ${extension})
    add_mex(ConstrainedRestartReadRestart ${type} ${include} ${extension})

    # Install all the script files
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
        DESTINATION share/optizelle/${type}/optizelle)

    # Copy in files for the unit tests.
    if(${unit})
        add_custom_target(
            ${type}_create_directory
            ALL
            COMMAND ${CMAKE_COMMAND} -E make_directory
                ${CMAKE_BINARY_DIR}/src/${type}/)
        add_custom_target(
            ${type}_mex_for_unit
            ALL
            COMMAND ${CMAKE_COMMAND} -E copy_directory
                ${CMAKE_BINARY_DIR}/src/matlab_octave/optizelle
                ${CMAKE_BINARY_DIR}/src/${type}/optizelle
            DEPENDS ${type}_create_directory)
        add_custom_target(
            ${type}_scripts_for_unit
            ALL
            COMMAND ${CMAKE_COMMAND} -E copy_directory
                ${CMAKE_SOURCE_DIR}/src/matlab_octave/optizelle
                ${CMAKE_BINARY_DIR}/src/${type}/optizelle
            DEPENDS ${type}_create_directory)
        add_custom_target(
            ${type}_setup_for_unit
            ALL
            COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_SOURCE_DIR}/src/matlab_octave/setupOptizelle.m
                ${CMAKE_BINARY_DIR}/src/${type}/
            DEPENDS ${type}_create_directory)
        add_custom_target(
            ${type}_jsonlab_for_unit
            ALL
            COMMAND ${CMAKE_COMMAND} -E copy_directory
                ${JSONLAB_DIR}
                ${CMAKE_BINARY_DIR}/src/thirdparty/share/jsonlab
            DEPENDS ${type}_create_directory)
    endif()
endmacro()

# Compile the MATLAB interface
if(ENABLE_MATLAB)
    # Necessary since we run through the unit test cmake scripts after this
    # directory where this variable is defined
    if(NOT ENABLE_MATLAB_UNIT)
        set(MATLAB_MEX_EXTENSION "MEX_UKNOWN")
        set(ENABLE_MATLAB_UNIT FALSE)
    endif()
    compile_interface(
        matlab
        ${MATLAB_INCLUDE_DIR}
        "${MATLAB_LIBRARY}"
        ${MATLAB_MEX_EXTENSION}
        ${ENABLE_MATLAB_UNIT})
endif()

# Compile the Octave interface
if(ENABLE_OCTAVE)
    if(NOT ENABLE_OCTAVE_UNIT)
        set(ENABLE_OCTAVE_UNIT FALSE)
    endif()
    compile_interface(
        octave
        ${OCTAVE_INCLUDE_DIR}
        "${OCTAVE_LIBRARY}"
        mex
        ${ENABLE_OCTAVE_UNIT})
endif()
