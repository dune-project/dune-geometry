install(FILES
  pointquadrature.hh
  simplexquadrature.hh
  genericquadrature.hh
  gauss_imp.hh
  jacobi_1_0_imp.hh
  jacobi_2_0_imp.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/geometry/quadraturerules)

#build the library libquadraturerules
dune_add_library(quadraturerules OBJECT
  cube.cc
  gauss.cc
  jacobi_1_0.cc
  jacobi_2_0.cc
  quadraturerules.cc
)


message(AUTHOR_WARNING "TODO in dune-geometry/dune/geometry/quadraturerules: Distribute jacobian.mac")
