install(FILES
  compositequadraturerule.hh
  pointquadrature.hh
  simplexquadrature.hh
  tensorproductquadrature.hh
  gaussquadrature.hh
  gausslobattoquadrature.hh
  gaussradauleftquadrature.hh
  gaussradaurightquadrature.hh
  jacobi1quadrature.hh
  jacobi2quadrature.hh
  jacobiNquadrature.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/geometry/quadraturerules)

exclude_from_headercheck(
  "pointquadrature.hh
  simplexquadrature.hh
  genericquadrature.hh
  gaussquadrature.hh
  gausslobattoquadrature.hh
  gaussradauleftquadrature.hh
  gaussradaurightquadrature.hh
  jacobi1quadrature.hh
  jacobi2quadrature.hh
  jacobiNquadrature.hh
  tensorproductquadrature.hh")

#build the library libquadraturerules
dune_add_library(quadraturerules OBJECT
  quadraturerules.cc
)
