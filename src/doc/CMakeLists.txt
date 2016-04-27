project(doc)

# Determine if we want an a4 sized document
mark_as_advanced(CLEAR ENABLE_A4_PAPER)
set(ENABLE_A4_PAPER OFF CACHE BOOL
    "Generate A4 sized documentation?  If false, the document is letter sized.")
if(ENABLE_A4_PAPER)
    set(PAPERTYPE "a4paper")
else()
    set(PAPERTYPE "letterpaper")
endif()

# Grab the LaTeX package
include(UseLATEX)

# Set a couple of paths to source files
set(OPTIZELLECPPPATH "${CMAKE_SOURCE_DIR}/src/cpp/optizelle/")
set(OPTIZELLEPYTHONPATH "${CMAKE_SOURCE_DIR}/src/python/Optizelle/")
set(OPTIZELLEMATLABPATH "${CMAKE_SOURCE_DIR}/src/matlab_octave/")
set(LICENSEPATH "${CMAKE_SOURCE_DIR}")
set(ROSENBROCKPATH "${CMAKE_SOURCE_DIR}/src/examples/rosenbrock/")
set(ROSENBROCKADVANCEDAPIPATH
    "${CMAKE_SOURCE_DIR}/src/examples/rosenbrock_advanced_api/")
set(SIMPLEEQUALITYPATH "${CMAKE_SOURCE_DIR}/src/examples/simple_equality/")
set(SIMPLEINEQUALITYPATH "${CMAKE_SOURCE_DIR}/src/examples/simple_inequality/")
set(SIMPLECONSTRAINEDPATH "${CMAKE_SOURCE_DIR}/src/examples/simple_constrained/")
set(SIMPLECONSTRAINEDADVANCEDAPIPATH
    "${CMAKE_SOURCE_DIR}/src/examples/simple_constrained_advanced_api/")
set(SIMPLEQUADRATICCONEPATH
    "${CMAKE_SOURCE_DIR}/src/examples/simple_quadratic_cone/")
set(SIMPLESDPCONEPATH
    "${CMAKE_SOURCE_DIR}/src/examples/simple_sdp_cone/")
set(RESTARTPATH "${CMAKE_SOURCE_DIR}/src/unit/restart/")
set(SQLRESTARTPATH "${CMAKE_SOURCE_DIR}/src/examples/sql_restart/")
set(CACHEPATH "${CMAKE_SOURCE_DIR}/src/examples/computation_caching/")
set(PROCESSEDPATH "${CMAKE_CURRENT_BINARY_DIR}")

# Compile the manual
add_latex_document(optizelle.tex DEFAULT_PDF
    CONFIGURE optizelle.tex
    USE_INDEX
    BIBFILES optizelle.bib
)

# Process the config file for documentation 
set(process_config_srcs "process_config.cpp" "file.cpp")
add_executable(process_config ${process_config_srcs})
add_custom_target(create_build_script
    COMMAND process_config
        "${CMAKE_SOURCE_DIR}/src/unit/config/config.cmake"
        "${PROCESSEDPATH}/mybuild.sh"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Turns a set of build flags into a build script"
    SOURCES ${process_config_srcs}
        "${CMAKE_SOURCE_DIR}/src/unit/config/config.cmake"
)
add_dependencies(pdf create_build_script)

# Run the Rosenbrock example and write the output to file
if(NOT ENABLE_CPP_UNIT AND NOT ENABLE_CPP_EXAMPLES)
    message(FATAL_ERROR "ENABLE_CPP_UNIT or ENABLE_CPP_EXAMPLES must be enabled in order to generate data for the manual")
endif()
add_custom_target(run_rosenbrock
    COMMAND rosenbrock
        "${CMAKE_SOURCE_DIR}/src/examples/rosenbrock/tr_newton.json"
        > "${PROCESSEDPATH}/rosen_output.txt"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Running the Rosenbrock example for the manual"
)
add_dependencies(pdf run_rosenbrock)

# Install the document with a nicer file name
install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/optizelle.pdf
    DESTINATION share/optizelle/doc
    RENAME Optizelle-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}.pdf)
