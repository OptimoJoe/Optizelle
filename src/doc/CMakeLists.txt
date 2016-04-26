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

# Compile the manual
add_latex_document(optizelle.tex DEFAULT_PDF
    CONFIGURE optizelle.tex
    USE_INDEX
    BIBFILES optizelle.bib
)
install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/optizelle.pdf
    DESTINATION share/optizelle/doc
    RENAME Optizelle-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}.pdf)
