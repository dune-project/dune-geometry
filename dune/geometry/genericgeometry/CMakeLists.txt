install(FILES
  codimtable.hh
  subtopologies.hh
  topologytypes.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/geometry/genericgeometry)

dune_add_library(genericgeometry OBJECT
  subtopologies.cc
  )
