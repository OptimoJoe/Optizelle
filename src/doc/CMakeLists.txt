project(doc)
add_latex_document(peopt.tex DEFAULT_PDF
    INPUTS rosen_code.tex rosen_input.tex rosen_output.tex
           simple_ineq_code.tex simple_ineq_input.tex simple_ineq_output.tex
           simple_eq_code.tex simple_eq_input.tex simple_eq_output.tex
           simple_con_code.tex simple_con_input.tex simple_con_output.tex
           objectiveSpecification.tex objectiveExample.tex
           constraintSpecification.tex constraintExample.tex
)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/peopt.pdf
    DESTINATION share/peopt/doc)
