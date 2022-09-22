# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

install(FILES
  base.cc
  hcube.cc
  hcubetriangulation.cc
  prismtriangulation.cc
  pyramidtriangulation.cc
  simplex.cc
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/geometry/refinement)
