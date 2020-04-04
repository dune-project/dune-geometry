set(HEADERS
  multilineargeometry.hh
  quadraturerules.hh
  referenceelements.hh
  type.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/python/geometry)
