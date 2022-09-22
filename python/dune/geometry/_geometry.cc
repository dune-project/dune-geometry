// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#include <dune/python/geometry/type.hh>
#include <dune/python/pybind11/pybind11.h>

PYBIND11_MODULE( _geometry, module )
{
  Dune::Python::registerGeometryType( module );

  // register geometry type contuctors
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
}
