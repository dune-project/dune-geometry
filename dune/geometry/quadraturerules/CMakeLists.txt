install(FILES
  calculatesimplexrules.hh
  compositequadraturerule.hh
  pointquadrature.hh
  simplexquadrature.hh
  tensorproductquadrature.hh
  gauss_imp.hh
  gausslobatto_imp.hh
  jacobi_1_0_imp.hh
  jacobi_2_0_imp.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/geometry/quadraturerules)

exclude_from_headercheck(
  "pointquadrature.hh
  simplexquadrature.hh
  genericquadrature.hh
  gauss_imp.hh
  gausslobatto_imp.hh
  jacobi_1_0_imp.hh
  jacobi_2_0_imp.hh
  tensorproductquadrature.hh
  calculatesimplexrules")

add_executable(intfunc intfunc.cc intfunc.hh calculatesimplexrules.hh)
#target_compile_definitions(intfunc PRIVATE -DUG_DIM_3)
target_compile_definitions(intfunc PRIVATE -DUG_DIM_3)
target_link_libraries(intfunc PRIVATE ${DUNE_LIBS})

#build the library libquadraturerules
dune_add_library(quadraturerules OBJECT
  gauss.cc
  jacobi_1_0.cc
  jacobi_2_0.cc
  quadraturerules.cc
  gausslobatto.cc
)
