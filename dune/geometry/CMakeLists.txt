add_subdirectory("quadraturerules")
add_subdirectory("genericgeometry")
add_subdirectory("refinement")
add_subdirectory("test"  EXCLUDE_FROM_ALL)

install(FILES
        affinegeometry.hh
        axisalignedcubegeometry.hh
        generalvertexorder.hh
        genericreferenceelements.hh
        mockgeometry.hh
	multilineargeometry.hh
        quadraturerules.hh
        referenceelements.hh
	refinement.hh
        topologyfactory.hh
        type.hh
        typeindex.hh
	virtualrefinement.hh
	virtualrefinement.cc
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/geometry
)

# install the header as done for the auto-tools
install(FILES test/checkgeometry.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/geometry/test)
