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
  foreach(DIM RANGE 0 4)
    add_library(registerrefelem_${DIM} OBJECT registerrefelem.cc)
    target_compile_definitions(registerrefelem_${DIM} PRIVATE DIM=${DIM})
    target_link_libraries(registerrefelem_${DIM} PRIVATE Dune::Geometry)

    # Register the object in the pybind11 module _geometry
    target_link_libraries(_geometry PRIVATE registerrefelem_${DIM})
  endforeach()
endif()

if(SKBUILD)
  install(TARGETS _geometry LIBRARY DESTINATION python/dune/geometry)
endif()
