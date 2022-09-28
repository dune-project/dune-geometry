# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

set(PICTURES gg_hexahedron.png gg_line.png gg_prism.png gg_pyramid.png
  gg_quadrilateral.png gg_tetrahedron.png gg_triangle.png
  gg_hexahedron_edges.png gg_tetrahedron_edges.png gg_prism_edges.png gg_pyramid_edges.png)

dune_create_inkscape_image_converter_target(
  TARGET referenceelement_images
  OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}
  IMAGES ${PICTURES}
)

# create reference element images before the documentation target
add_dependencies(doc referenceelement_images)

install(FILES ${PICTURES} DESTINATION ${CMAKE_INSTALL_DOCDIR}/appl/refelements)
