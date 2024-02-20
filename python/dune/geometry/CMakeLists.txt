# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

add_python_targets(geometry
  __init__
  _referenceelements
  quadpy
)

set(DIMS 0 1 2 3 4)

set(GEOMETRY_SOURCES _geometry.cc)
if( DUNE_ENABLE_PYTHONMODULE_PRECOMPILE )
  foreach( DIM ${DIMS} )
    set(GEOMETRY_SOURCES ${GEOMETRY_SOURCES} registerrefelem_${DIM}.cc)
  endforeach()
endif()

dune_add_pybind11_module(NAME _geometry SOURCES ${GEOMETRY_SOURCES} )
set_property(TARGET _geometry PROPERTY LINK_LIBRARIES Dune::Geometry APPEND)

if(SKBUILD)
  install(TARGETS _geometry LIBRARY DESTINATION python/dune/geometry)
endif()
