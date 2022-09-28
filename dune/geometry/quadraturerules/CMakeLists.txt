# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

install(FILES
  compositequadraturerule.hh
  gausslobattoquadrature.hh
  gaussquadrature.hh
  gaussradauleftquadrature.hh
  gaussradaurightquadrature.hh
  jacobi1quadrature.hh
  jacobi2quadrature.hh
  jacobiNquadrature.hh
  pointquadrature.hh
  prismquadrature.hh
  simplexquadrature.hh
  tensorproductquadrature.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/geometry/quadraturerules)

exclude_from_headercheck(
  "gausslobattoquadrature.hh
  gaussquadrature.hh
  gaussradauleftquadrature.hh
  gaussradaurightquadrature.hh
  jacobi1quadrature.hh
  jacobi2quadrature.hh
  jacobiNquadrature.hh
  pointquadrature.hh
  prismquadrature.hh
  simplexquadrature.hh
  tensorproductquadrature.hh")

#build the library libquadraturerules
target_sources(dunegeometry PRIVATE
  quadraturerules.cc
)
