set(TESTS
  test-affinegeometry
  test-axisalignedcubegeometry
  test-geometrytype
  test-multilineargeometry
  test-quadrature
  test-referenceelements
  test-refinement)

add_executable(test-affinegeometry test-affinegeometry.cc)
add_executable(test-axisalignedcubegeometry test-axisalignedcubegeometry.cc)
add_executable(test-geometrytype test-geometrytype.cc)
add_executable(test-referenceelements test-referenceelements.cc)
add_executable(test-quadrature test-quadrature.cc)
add_executable(test-multilineargeometry test-multilineargeometry.cc)
add_executable(test-refinement test-refinement.cc)

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS})

foreach(_test ${TESTS})
  target_link_libraries(${_test} dunegeometry ${DUNE_LIBS})
  add_test(${_test} ${_test})
endforeach(_test ${TESTS})
