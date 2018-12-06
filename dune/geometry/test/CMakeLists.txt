dune_add_test(SOURCES test-affinegeometry.cc
              LINK_LIBRARIES dunegeometry)

dune_add_test(SOURCES test-axisalignedcubegeometry.cc
              LINK_LIBRARIES dunegeometry)

# whether the compiler is affected by
# https://gitlab.dune-project.org/core/dune-geometry/issues/19,
# a.k.a. https://gcc.gnu.org/bugzilla/show_bug.cgi?id=87288
set(ISSUE19_UNAFFECTED_GXX ON)
if((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") AND
    (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 8.2) AND
    (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9))
  set(ISSUE19_UNAFFECTED_GXX OFF)
endif()
dune_add_test(SOURCES test-cornerstoragerefwrap.cc
              LINK_LIBRARIES dunegeometry
              CMAKE_GUARD ISSUE19_UNAFFECTED_GXX)

dune_add_test(SOURCES test-fromvertexcount.cc)

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

dune_add_test(SOURCES test-constexpr-geometrytype.cc
              LINK_LIBRARIES dunegeometry)
