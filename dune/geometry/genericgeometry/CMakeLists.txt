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

# Install some test headers, because they get used by tests in other modules
# We do this here as test will not be considered for make install
install(FILES test/checkgeometry.hh DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/geometry/genericgeometry/test)
