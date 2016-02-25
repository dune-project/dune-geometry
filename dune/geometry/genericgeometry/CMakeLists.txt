add_subdirectory("test")
install(FILES
  codimtable.hh
  conversion.hh
  geometrytraits.hh
  matrixhelper.hh
  referencedomain.hh
  referenceelements.hh
  subtopologies.hh
  topologytypes.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/geometry/genericgeometry)

dune_add_library(genericgeometry OBJECT
  referencedomain.cc
  subtopologies.cc
  )
