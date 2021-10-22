add_python_targets(geometry
  __init__
  _referenceelements
  quadpy
)
dune_add_pybind11_module(NAME _geometry)
set_property(TARGET _geometry PROPERTY LINK_LIBRARIES dunecommon dunegeometry APPEND)
if(SKBUILD)
  install(TARGETS _geometry LIBRARY DESTINATION python/dune/geometry)
endif()
