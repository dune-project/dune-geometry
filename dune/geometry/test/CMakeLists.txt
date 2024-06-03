# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

dune_add_test(SOURCES benchmark-geometries.cc
              LINK_LIBRARIES dunegeometry)

dune_add_test(SOURCES test-affinegeometry.cc
              LINK_LIBRARIES dunegeometry)

dune_add_test(SOURCES test-axisalignedcubegeometry.cc
              LINK_LIBRARIES dunegeometry)

dune_add_test(SOURCES test-cornerstoragerefwrap.cc
              LINK_LIBRARIES dunegeometry)

dune_add_test(SOURCES test-fromvertexcount.cc)

dune_add_test(SOURCES test-referenceelements.cc
              LINK_LIBRARIES dunegeometry)

dune_add_test(SOURCES test-quadrature.cc
              LINK_LIBRARIES dunegeometry)

dune_add_test(SOURCES test-mappedgeometry.cc
              LINK_LIBRARIES dunegeometry)

dune_add_test(SOURCES test-multilineargeometry.cc
              LINK_LIBRARIES dunegeometry)

dune_add_test(SOURCES test-nonetype.cc
              LINK_LIBRARIES dunegeometry)

dune_add_test(SOURCES test-localfiniteelementgeometry.cc
              LINK_LIBRARIES dunegeometry)

dune_add_test(SOURCES test-refinement.cc
              LINK_LIBRARIES dunegeometry)

dune_add_test(SOURCES test-constexpr-geometrytype.cc
              LINK_LIBRARIES dunegeometry)

dune_add_test(SOURCES test-geometrytype-id.cc
              LINK_LIBRARIES dunegeometry)
