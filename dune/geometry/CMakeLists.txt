# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

add_subdirectory("quadraturerules")
add_subdirectory("refinement")
add_subdirectory("utility")
add_subdirectory("test")

install(FILES
  affinegeometry.hh
  axisalignedcubegeometry.hh
  dimension.hh
  generalvertexorder.hh
  multilineargeometry.hh
  quadraturerules.hh
  referenceelement.hh
  referenceelementimplementation.hh
  referenceelements.hh
  refinement.hh
  topologyfactory.hh
  type.hh
  typeindex.hh
  virtualrefinement.hh
  virtualrefinement.cc
  deprecated_topology.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/geometry
)

# install the header as done for the auto-tools
install(FILES test/checkgeometry.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/geometry/test)

target_sources(dunegeometry PRIVATE
  referenceelementimplementation.cc
  )
