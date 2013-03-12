
if(INKSCAPE)
  set(DOCS ${DOCS} referencetokuhn0.png simplexvertexindex.png)
  set(EPS_IMGS simplexvertexindex.eps referencetokuhn0.eps kuhntriangulation.eps)
  inkscape_generate_eps_from_svg(DPI 180
    OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}
    INPUT_DIR ${CMAKE_CURRENT_BINARY_DIR} ${EPS_IMGS})
endif(INKSCAPE)
dune_add_latex_document(refinement.tex FATHER_TARGET doc
  DEFAULT_SAFEPDF
  DEPENDS ${EPS_IMGS})
set(PDFDOCS ${CMAKE_CURRENT_BINARY_DIR}/refinement.pdf)

create_doc_install(${CMAKE_CURRENT_BINARY_DIR}/refinement.pdf ${CMAKE_INSTALL_DOCDIR}/refinement)
