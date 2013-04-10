set(TESTS
  testsubtopologies
  testbasicgeometry)

set(GEOMETRYTYPE "cube")
set(DIMENSION "3")

add_executable(testsubtopologies testsubtopologies.cc)
add_executable(testbasicgeometry testbasicgeometry.cc)
target_link_libraries(testbasicgeometry dunecommon dunegeometry)
add_executable(testmatrix testmatrix.cc)

set_property(TARGET testsubtopologies APPEND_STRING  PROPERTY COMPILE_FLAGS
" -DGEOMETRYTYPE=${GEOMETRYTYPE} -DDIMENSION=${DIMENSION}")
add_dune_gmp_flags(testmatrix)

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)

add_dependencies(${_test_target} ${TESTS})
foreach(_test ${TESTS})
  add_test(${_test} ${_test})
endforeach(_test ${TESTS})
