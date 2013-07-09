add_subdirectory("test" EXCLUDE_FROM_ALL)
install(FILES
        cachedmapping.hh
	codimtable.hh
	conversion.hh
        cornermapping.hh
	geometrytraits.hh
        geometry.hh
	hybridmapping.hh
        mapping.hh
	mappingprovider.hh
	matrixhelper.hh
        maximum.hh
	referencedomain.hh
	referenceelements.hh
        subtopologies.hh
        topologytypes.hh
	traceprovider.hh
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/geometry/genericgeometry)

dune_add_library(genericgeometry OBJECT
  referencedomain.cc
  subtopologies.cc
  )
