install(FILES
  compositequadraturerule.hh
  nocopyvector.hh
  pointquadrature.hh
  gaussquadrature.inc.hh
  gausslobattoquadrature.inc.hh
  jacobi1quadrature.inc.hh
  jacobi2quadrature.inc.hh
  cube2dquadrature.inc.hh
  cube3dquadrature.inc.hh
  simplex2dquadrature.inc.hh
  simplex3dquadrature.inc.hh
  prismquadrature.inc.hh
  pyramidquadrature.inc.hh
  tensorproductquadrature.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/geometry/quadraturerules)

exclude_from_headercheck(
  "pointquadrature.hh
  genericquadrature.hh
  gaussquadrature.inc.hh
  gausslobattoquadrature.inc.hh
  jacobi1quadrature.inc.hh
  jacobi2quadrature.inc.hh
  cube2dquadrature.inc.hh
  cube3dquadrature.inc.hh
  simplex2dquadrature.inc.hh
  simplex3dquadrature.inc.hh
  prismquadrature.inc.hh
  pyramidquadrature.inc.hh
  tensorproductquadrature.hh")

#build the library libquadraturerules
dune_add_library(quadraturerules OBJECT
  quadraturerules.cc
)
