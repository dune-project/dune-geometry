install(FILES
  codimtable.hh
  referencedomain.hh
  subtopologies.hh
  topologytypes.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/geometry/genericgeometry)

dune_add_library(genericgeometry OBJECT
  referencedomain.cc
  subtopologies.cc
  )
