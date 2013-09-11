project(doc)
include(UseLATEX)
add_latex_document(optizelle.tex DEFAULT_PDF
    INPUTS rosen_code.tex rosen_input.tex rosen_output.tex
           simple_ineq_code.tex simple_ineq_input.tex simple_ineq_output.tex
           simple_eq_code.tex simple_eq_input.tex simple_eq_output.tex
           simple_con_code.tex simple_con_input.tex simple_con_output.tex
           objectiveSpecification.tex objectiveExample.tex
           constraintSpecification.tex constraintExample.tex
           initialGuessExample.tex
           allocateCodomainExample.tex
           stateSpecification.tex
           jsonReaderSpecification.tex jsonReaderExample1.tex
               jsonReaderExample2.tex
           functionBundleSpecification.tex
           solverSpecification.tex
           solutionExtractionSpecification.tex
           messagingSpecification.tex
           vectorSpaceSpecification.tex
)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/optizelle.pdf
    DESTINATION share/optizelle/doc)
