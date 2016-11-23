dune_add_test(SOURCES test-affinegeometry.cc
              LINK_LIBRARIES dunegeometry)

dune_add_test(SOURCES test-axisalignedcubegeometry.cc
              LINK_LIBRARIES dunegeometry)

dune_add_test(SOURCES test-cornerstoragerefwrap.cc
              LINK_LIBRARIES dunegeometry)

dune_add_test(SOURCES test-geometrytype.cc)

dune_add_test(SOURCES test-referenceelements.cc
              LINK_LIBRARIES dunegeometry)

dune_add_test(SOURCES test-quadrature.cc
              LINK_LIBRARIES dunegeometry)

dune_add_test(SOURCES test-multilineargeometry.cc
              LINK_LIBRARIES dunegeometry)

dune_add_test(SOURCES test-nonetype.cc
              LINK_LIBRARIES dunegeometry)

dune_add_test(SOURCES test-refinement.cc
              LINK_LIBRARIES dunegeometry)
