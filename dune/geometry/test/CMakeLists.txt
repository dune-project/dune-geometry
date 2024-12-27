# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

# Link all test targets in this directory against Dune::Geometry
link_libraries(Dune::Geometry)

dune_add_test(SOURCES benchmark-geometries.cc)

dune_add_test(SOURCES test-affinegeometry.cc)

dune_add_test(SOURCES test-axisalignedcubegeometry.cc)

dune_add_test(SOURCES test-cornerstoragerefwrap.cc)

dune_add_test(SOURCES test-fromvertexcount.cc)

dune_add_test(SOURCES test-referenceelements.cc)

dune_add_test(SOURCES test-quadrature.cc)

dune_add_test(SOURCES test-mappedgeometry.cc)

dune_add_test(SOURCES test-multilineargeometry.cc)

dune_add_test(SOURCES test-nonetype.cc)

dune_add_test(SOURCES test-localfiniteelementgeometry.cc)

dune_add_test(SOURCES test-refinement.cc)

dune_add_test(SOURCES test-constexpr-geometrytype.cc)

dune_add_test(SOURCES test-geometrytype-id.cc)
