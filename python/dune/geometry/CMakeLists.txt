# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

add_python_targets(geometry
  __init__
  _referenceelements
  quadpy
)


dune_add_pybind11_module(NAME _geometry SOURCES _geometry.cc)
target_link_libraries(_geometry PUBLIC Dune::Geometry)
if( DUNE_ENABLE_PYTHONMODULE_PRECOMPILE )
  set(DIMS 0 1 2 3 4)
  foreach( DIM ${DIMS} )
    dune_add_pybind11_submodule(NAME registerrefelem_${DIM}
                                MODULE _geometry SOURCES registerrefelem.cc
                                COMPILE_DEFINITIONS DIM=${DIM})
  endforeach()
endif()
set_property(TARGET _geometry PROPERTY LINK_LIBRARIES Dune::Geometry APPEND)

if(SKBUILD)
  install(TARGETS _geometry LIBRARY DESTINATION python/dune/geometry)
endif()
