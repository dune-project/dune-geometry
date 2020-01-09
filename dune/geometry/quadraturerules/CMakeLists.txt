install(FILES
  compositequadraturerule.hh
  pointquadrature.hh
  simplexquadrature.hh
  tensorproductquadrature.hh
  gauss_imp.hh
  gausslobatto_imp.hh
  jacobi_1_0_imp.hh
  jacobi_2_0_imp.hh
  jacobi_n_0.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/geometry/quadraturerules)

exclude_from_headercheck(
  "pointquadrature.hh
  simplexquadrature.hh
  genericquadrature.hh
  gauss_imp.hh
  gausslobatto_imp.hh
  jacobi_1_0_imp.hh
  jacobi_2_0_imp.hh
  jacobi_n_0.hh
  tensorproductquadrature.hh")

#build the library libquadraturerules
dune_add_library(quadraturerules OBJECT
  gauss.cc
  jacobi_1_0.cc
  jacobi_2_0.cc
  quadraturerules.cc
  gausslobatto.cc
)
