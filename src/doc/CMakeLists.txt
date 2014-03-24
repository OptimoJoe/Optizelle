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
set(LICENSEPATH "${CMAKE_SOURCE_DIR}")
set(ROSENBROCKPATH "${CMAKE_SOURCE_DIR}/src/examples/rosenbrock/")
set(SIMPLEEQUALITYPATH "${CMAKE_SOURCE_DIR}/src/examples/simple_equality/")
set(OUTPUTPATH "${CMAKE_CURRENT_SOURCE_DIR}")
set(OPTIZELLECPPPATH "${CMAKE_SOURCE_DIR}/src/cpp/optizelle/")
set(OPTIZELLEPYTHONPATH "${CMAKE_SOURCE_DIR}/src/python/Optizelle/")
set(OPTIZELLEMATLABPATH "${CMAKE_SOURCE_DIR}/src/matlab/")

# Compile the manual
add_latex_document(optizelle.tex DEFAULT_PDF
    INPUTS messagingSpecification.tex
           vectorSpaceSpecification.tex
    CONFIGURE optizelle.tex
    USE_INDEX
)
install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/optizelle.pdf
    DESTINATION share/optizelle/doc
    RENAME Optizelle-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}.pdf)
