install(FILES
  cube_imp.hh
  gaussquadrature.hh
  genericquadrature.hh
  jacobi_1_0_imp.hh
  jacobi_2_0_imp.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIRS}/dune/geometry/quadraturerules)

#build the library libquadraturerules
dune_add_library("quadraturerules"
  cube.cc
  jacobi_1_0.cc
  jacobi_2_0.cc
  quadraturerules.cc
)
dune_target_link_libraries("quadraturerules" dunecommon)

message(AUTHOR_WARNING "TODO in dune-geometry/dune/geometry/quadraturerules: Distribute jacobian.mac")
