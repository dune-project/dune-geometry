add_subdirectory(dune)

dune_python_install_package(
  PATH "."
  CMAKE_METADATA_FILE dune/geometry/metadata.cmake
  DEPENDS _geometry
)
