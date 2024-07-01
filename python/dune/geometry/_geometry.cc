// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#include <string>

#include <dune-common-config.hh> // DUNE_ENABLE_PYTHONMODULE_PRECOMPILE
#include <dune/python/geometry/type.hh>
#include <dune/python/geometry/referenceelements.hh>
#include <dune/python/pybind11/pybind11.h>

#ifdef DUNE_ENABLE_PYTHONMODULE_PRECOMPILE
#include "registerrefelem.hh"
#endif

PYBIND11_MODULE( _geometry, module )
{
  Dune::Python::registerGeometryType( module );

  // register geometry type constructors
  module.def( "simplex", [] ( int dim ) { return Dune::GeometryTypes::simplex( dim ); } );
  module.def( "cube", [] ( int dim ) { return Dune::GeometryTypes::cube( dim ); } );
  module.def( "none", [] ( int dim ) { return Dune::GeometryTypes::none( dim ); } );

  // register predefined geometry types
  module.attr( "vertex" ) = Dune::GeometryTypes::vertex;
  module.attr( "line" ) = Dune::GeometryTypes::line;
  module.attr( "triangle" ) = Dune::GeometryTypes::triangle;
  module.attr( "quadrilateral" ) = Dune::GeometryTypes::quadrilateral;
  module.attr( "tetrahedron" ) = Dune::GeometryTypes::tetrahedron;
  module.attr( "pyramid" ) = Dune::GeometryTypes::pyramid;
  module.attr( "prism" ) = Dune::GeometryTypes::prism;
  module.attr( "hexahedron" ) = Dune::GeometryTypes::hexahedron;

#ifdef DUNE_ENABLE_PYTHONMODULE_PRECOMPILE
  registerReferenceElementSubModule<0>(module);
  registerReferenceElementSubModule<1>(module);
  registerReferenceElementSubModule<2>(module);
  registerReferenceElementSubModule<3>(module);
  registerReferenceElementSubModule<4>(module);
#endif
}
