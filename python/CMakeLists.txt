add_subdirectory(dune)

dune_python_install_package(
  PATH "."
  CMAKE_METADATA_FILE dune/geometry/metadata.cmake
  DEPENDS _geometry
  CMAKE_METADATA_FLAGS DUNE_OPTS_FILE
)
if(POLICY CMP0087)
  cmake_policy(SET CMP0087 NEW)
endif()
