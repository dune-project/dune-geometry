# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

dune_add_latex_document(
  SOURCE refinement.tex
  TARGET dune_geometry_refinement_tex
  FATHER_TARGET doc
  INSTALL ${CMAKE_INSTALL_DOCDIR})

# dune_add_latex_document doesn't create a target if latex is not found
if(TARGET dune_geometry_refinement_tex)
  set(PNG_IMGS simplexvertexindex.png referencetokuhn0.png kuhntriangulation.png)
  dune_create_inkscape_image_converter_target(
    TARGET refinement_images
    DPI 180
    IMAGES ${PNG_IMGS}
  )

  # create images before latex document
  add_dependencies(dune_geometry_refinement_tex refinement_images)
endif()
