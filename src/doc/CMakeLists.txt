project(doc)
add_latex_document(peopt.tex DEFAULT_PDF)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/peopt.pdf
    DESTINATION share/peopt/doc)
