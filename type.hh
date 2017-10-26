// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PYTHON_GEOMETRY_TYPE_HH
#define DUNE_PYTHON_GEOMETRY_TYPE_HH

#include <dune/geometry/type.hh>

#include <dune/python/pybind11/operators.h>
#include <dune/python/pybind11/pybind11.h>

namespace Dune
{

  namespace Python
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

      cls.def( "__str__", [] ( const GeometryType &self ) -> std::string {
          if( self.isNone() )
            return "none(" + std::to_string( self.dim() ) + ")";
          switch( self.dim() )
          {
          case 0:
            return "vertex";
          case 1:
            return "line";
          case 2:
            return (self.isSimplex() ? "triangle" : "quadrilateral");
          case 3:
            if( self.isSimplex() )
              return "tetrahedron";
            else if( self.isHexahedron() )
              return "hexahedron";
            else if( self.isPyramid() )
              return "pyramid";
            else if( self.isPrism() )
              return "prism";
          default:
            if( self.isSimplex() )
              return "simplex(" + std::to_string( self.dim() ) + ")";
            else if( self.isCube() )
              return "cube(" + std::to_string( self.dim() ) + ")";
            else
              return "general(" + std::to_string( self.id() ) + ", " + std::to_string( self.dim() ) + ")";
          }
        } );

      cls.def_property_readonly( "dim", [] ( const GeometryType &self ) { return self.dim(); } );
      cls.def_property_readonly( "dimension", [] ( const GeometryType &self ) { return self.dim(); } );

      return cls;
    }

  } // namespace Python

} // namespace Dune

#endif // ifndef DUNE_PYTHON_GEOMETRY_TYPE_HH
