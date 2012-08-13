set(PICTURES hexahedron.png line.png prism.png pyramid.png
  quadrilateral.png tetrahedron.png triangle.png
  hexahedron_edges.png tetrahedron_edges.png prism_edges.png pyramid_edges.png)

inkscape_generate_png_from_svg(${PICTURES} OUTPUT_DIR ${CMAKE_BINARY_DIR})
add_custom_target(update_images DEPENDS ${PICTURES})

install(FILES ${PICTURES} DESTINATION ${CMAKE_INSTALL_DOCDIR}/appl/refelements)
