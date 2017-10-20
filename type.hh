// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PYTHON_GEOMETRY_TYPE_HH
#define DUNE_PYTHON_GEOMETRY_TYPE_HH

#include <dune/geometry/type.hh>

#include <dune/python/pybind11/operators.h>
#include <dune/python/pybind11/pybind11.h>

namespace Dune
{

  namespace CorePy
  {

    pybind11::class_< GeometryType > registerGeometryType( pybind11::module scope )
    {
      pybind11::class_< GeometryType > cls( scope, "GeometryType" );

      cls.def_property_readonly( "isVertex", [] ( const GeometryType &self ) { return self.isVertex(); } );
      cls.def_property_readonly( "isLine", [] ( const GeometryType &self ) { return self.isLine(); } );
      cls.def_property_readonly( "isTriangle", [] ( const GeometryType &self ) { return self.isTriangle(); } );
      cls.def_property_readonly( "isQuadrilateral", [] ( const GeometryType &self ) { return self.isQuadrilateral(); } );
      cls.def_property_readonly( "isTetrahedron",[] ( const GeometryType &self ) { return self.isTetrahedron(); } );
      cls.def_property_readonly( "isPyramid",[] ( const GeometryType &self ) { return self.isPyramid(); } );
      cls.def_property_readonly( "isPrism", [] ( const GeometryType &self ) { return self.isPrism(); } );
      cls.def_property_readonly( "isHexahedron", [] ( const GeometryType &self ) { return self.isHexahedron(); } );
      cls.def_property_readonly( "isSimplex", [] ( const GeometryType &self ) { return self.isSimplex(); } );
      cls.def_property_readonly( "isCube", [] ( const GeometryType &self ) { return self.isCube(); } );
      cls.def_property_readonly( "isNone", [] ( const GeometryType &self ) { return self.isNone(); } );

      cls.def( pybind11::self == pybind11::self );
      cls.def( pybind11::self != pybind11::self );
      cls.def( "__hash__", [] ( const GeometryType &self ) { return self.id(); } );

      cls.def_property_readonly( "dim", [] ( const GeometryType &self ) { return self.dim(); } );
      cls.def_property_readonly( "dimension", [] ( const GeometryType &self ) { return self.dim(); } );

      return cls;
    }

  } // namespace CorePy

} // namespace Dune

#endif // ifndef DUNE_PYTHON_GEOMETRY_TYPE_HH
