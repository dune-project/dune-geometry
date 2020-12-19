add_python_targets(geometry
  __init__
  _referenceelements
  quadpy
)

dune_add_pybind11_module(NAME _geometry)
set_property(TARGET _geometry PROPERTY LINK_LIBRARIES dunecommon dunegeometry APPEND)
install(TARGETS _geometry LIBRARY DESTINATION python/dune/geometry)

set( DIMS 0 1 2 3 )
foreach( DIM ${DIMS} )
  dune_add_pybind11_module(NAME _defaultreferenceelements_${DIM})
  set_property(TARGET _defaultreferenceelements_${DIM} PROPERTY LINK_LIBRARIES dunecommon dunegeometry APPEND)
  install(TARGETS _defaultreferenceelements_${DIM} LIBRARY DESTINATION python/dune/geometry)
endforeach()
