set(TESTS
  test_affinegeometry
  test_axisalignedcubegeometry
  test_geometrytype
  test_multilineargeometry
  test_quadrature
  test_referenceelements
  test_refinement)

add_executable(test_affinegeometry test-affinegeometry.cc)
target_link_libraries(test_affinegeometry dunecommon dunegeometry)
add_executable(test_axisalignedcubegeometry test-axisalignedcubegeometry.cc)
target_link_libraries(test_axisalignedcubegeometry dunecommon dunegeometry)
add_executable(test_geometrytype "test-geometrytype.cc")
target_link_libraries(test_geometrytype dunecommon dunegeometry)
add_executable(test_referenceelements "test-referenceelements.cc")
target_link_libraries(test_referenceelements dunegeometry)
add_executable(test_quadrature "test-quadrature.cc")
target_link_libraries(test_quadrature dunecommon dunegeometry)
add_executable(test_multilineargeometry test-multilineargeometry.cc)
target_link_libraries(test_multilineargeometry dunecommon dunegeometry)
add_executable(test_refinement test-refinement.cc)
target_link_libraries(test_refinement dunecommon dunegeometry)

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS})

foreach(_test ${TESTS})
  add_test(${_test} ${_test})
endforeach(_test ${TESTS})
