set(TESTS
  testsubtopologies
  testbasicgeometry)

set(GEOMETRYTYPE "cube")
set(DIMENSION "3")

add_executable(testsubtopologies testsubtopologies.cc)
add_executable(testbasicgeometry testbasicgeometry.cc)
add_executable(testmatrix testmatrix.cc)

set_property(TARGET testsubtopologies APPEND_STRING  PROPERTY COMPILE_FLAGS
" -DGEOMETRYTYPE=${GEOMETRYTYPE} -DDIMENSION=${DIMENSION}")

message(AUTHOR_WARNING "We have not added GMP flags in dune/geomerty/genericgeomerty/test as there is no test activated for it in the trunk")
