set(PICTURES gg_hexahedron.png gg_line.png gg_prism.png gg_pyramid.png
  gg_quadrilateral.png gg_tetrahedron.png gg_triangle.png
  gg_hexahedron_edges.png gg_tetrahedron_edges.png gg_prism_edges.png gg_pyramid_edges.png)

inkscape_generate_png_from_svg(${PICTURES} OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR})
add_custom_target(update_images DEPENDS ${PICTURES})

install(FILES ${PICTURES} DESTINATION ${CMAKE_INSTALL_DOCDIR}/appl/refelements)
