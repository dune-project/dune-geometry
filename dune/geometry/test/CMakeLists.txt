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
    (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8.3))
  separate_arguments(opts UNIX_COMMAND "${CMAKE_CXX_FLAGS}")
  execute_process(
    COMMAND ${CMAKE_CXX_COMPILER} ${opts} -Q --help=optimizers
    OUTPUT_VARIABLE opt_settings
    ERROR_VARIABLE ignored # silence cc1 complaint about c++ options...
  )
  if(opt_settings MATCHES "-ftree-loop-vectorize[ \t]+\\[enabled\\]")
    message(WARNING "You are using gcc 8.2 with optimization options "
      "including -ftree-loop-vectorize (possibly implicated by -O3).  This is "
      "known to miscompile, see "
      "https://gitlab.dune-project.org/core/dune-geometry/issues/19 and "
      "https://gcc.gnu.org/bugzilla/show_bug.cgi?id=87288")
    message(WARNING "Skipping test-cornerstoragerefwrap due to known-bad "
      "compiler optimization combination")
    set(ISSUE19_UNAFFECTED_GXX OFF)
  endif()
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
