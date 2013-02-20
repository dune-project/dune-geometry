add_subdirectory("quadraturerules")
add_subdirectory("genericgeometry")
add_subdirectory("test"  EXCLUDE_FROM_ALL)

install(FILES
        generalvertexorder.hh
        genericreferenceelements.hh
        mockgeometry.hh
        quadraturerules.hh
        referenceelements.hh
        topologyfactory.hh
        type.hh
        typeindex.hh
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/geometry
)
