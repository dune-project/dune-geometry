if(INKSCAPE)
  set(EPS_IMGS simplexvertexindex.eps referencetokuhn0.eps kuhntriangulation.eps)
  inkscape_generate_eps_from_svg(DPI 180
    OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}
    INPUT_DIR ${CMAKE_CURRENT_BINARY_DIR} ${EPS_IMGS})
endif()

dune_add_latex_document(refinement.tex
  DEPENDS ${EPS_IMGS})

create_doc_install(${CMAKE_CURRENT_BINARY_DIR}/refinement.pdf ${CMAKE_INSTALL_DOCDIR}/refinement refinement_safepdf)
