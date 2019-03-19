set(PNG_IMGS simplexvertexindex.png referencetokuhn0.png kuhntriangulation.png)
inkscape_generate_png_from_svg(
 DPI 180
 ${PNG_IMGS})

dune_add_latex_document(
  SOURCE refinement.tex
  FATHER_TARGET doc
  INSTALL ${CMAKE_INSTALL_DOCDIR})
