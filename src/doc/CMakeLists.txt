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

# Compile the manual.  Note, this file is basically just a stub that includes
# the prepocessed and configured optizelle.tex below.
add_latex_document(master.tex
    CONFIGURE master.tex
    DEFAULT_PDF
    USE_INDEX
    BIBFILES optizelle.bib
)

# Process the config file for documentation 
preprocess(optizelle.tex "${CMAKE_SOURCE_DIR}/src/unit/config/flags.txt")
add_custom_target(configure_paths
    COMMAND ${CMAKE_COMMAND}
        -D MY_SOURCE_DIR:PATH=${CMAKE_SOURCE_DIR}
        -D MY_BINARY_DIR:PATH=${CMAKE_SOURCE_DIR}
        -D MY_CURRENT_BINARY_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}
        -D CPACK_PACKAGE_VERSION_MAJOR=${CPACK_PACKAGE_VERSION_MAJOR}
        -D CPACK_PACKAGE_VERSION_MINOR=${CPACK_PACKAGE_VERSION_MINOR}
        -D CPACK_PACKAGE_VERSION_PATCH=${CPACK_PACKAGE_VERSION_PATCH}
        -P ${CMAKE_CURRENT_SOURCE_DIR}/CMakePathConfigure.txt)
add_dependencies(configure_paths optizelle.tex.pp.run)
add_dependencies(pdf configure_paths)

# Run the Rosenbrock example and write the output to file
if(NOT ENABLE_CPP_UNIT AND NOT ENABLE_CPP_EXAMPLES)
    message(FATAL_ERROR "ENABLE_CPP_UNIT or ENABLE_CPP_EXAMPLES must be enabled in order to generate data for the manual")
endif()
add_custom_target(run_rosenbrock
    COMMAND rosenbrock
        "${CMAKE_SOURCE_DIR}/src/examples/rosenbrock/tr_newton.json"
        > "${CMAKE_CURRENT_BINARY_DIR}/rosen_output.txt"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Running the Rosenbrock example for the manual"
)
add_dependencies(pdf run_rosenbrock)

# Install the document with a nicer file name
install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/master.pdf
    DESTINATION share/optizelle/doc
    RENAME Optizelle-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}.pdf)
