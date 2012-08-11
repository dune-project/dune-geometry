set(TEST
  test-quadrature
  test-genericquadrature
  test-referenceelements)

add_executable(test-referenceelements "test-referenceelements.cc")
add_executable(test-quadrature "test-quadrature.cc")
add_executable(test-genericquadrature "test-genericquadrature.cc")
