set(TESTS
  test-axisalignedcubegeometry
  test-genericquadrature
  test-quadrature
  test-referenceelements)

add_executable(test-axisalignedcubegeometry test-axisalignedcubegeometry.cc)
target_link_libraries(test-axisalignedcubegeometry dunecommon dunegeometry)
add_executable(test-referenceelements "test-referenceelements.cc")
add_executable(test-quadrature "test-quadrature.cc")
target_link_libraries(test-quadrature dunecommon dunegeometry)
add_executable(test-genericquadrature "test-genericquadrature.cc")
target_link_libraries(test-genericquadrature dunecommon dunegeometry)

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS})

foreach(_test ${TESTS})
  add_test(${_test} ${_test})
endforeach(_test ${TESTS})
