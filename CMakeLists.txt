add_python_targets(geometry
  __init__
)

add_library(_geometry SHARED _geometry.cc)
set_target_properties(_geometry PROPERTIES PREFIX "")
