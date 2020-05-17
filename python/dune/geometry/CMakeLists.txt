add_python_targets(geometry
  __init__
  _referenceelements
  _geometry
  quadpy
)
dune_symlink_to_source_files(FILES _geometry.cc)
