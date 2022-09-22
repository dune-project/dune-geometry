// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_PYTHON_GEOMETRY_TYPE_HH
#define DUNE_PYTHON_GEOMETRY_TYPE_HH

#include <cassert>

#include <dune/common/exceptions.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

#include <dune/python/pybind11/operators.h>
#include <dune/python/pybind11/pybind11.h>

namespace Dune
{

  // to_string for GeometryType
  // --------------------------

  inline static std::string to_string ( const GeometryType &type )
  {
    if( type.isNone() )
      return "none(" + std::to_string( type.dim() ) + ")";
    switch( type.dim() )
    {
    case 0:
      return "vertex";
    case 1:
      return "line";
    case 2:
      return (type.isSimplex() ? "triangle" : "quadrilateral");
    case 3:
      if( type.isSimplex() )
        return "tetrahedron";
      else if( type.isHexahedron() )
        return "hexahedron";
      else if( type.isPyramid() )
        return "pyramid";
      else if( type.isPrism() )
        return "prism";
    default:
      if( type.isSimplex() )
        return "simplex(" + std::to_string( type.dim() ) + ")";
      else if( type.isCube() )
        return "cube(" + std::to_string( type.dim() ) + ")";
      else
        return "general(" + std::to_string( type.id() ) + ", " + std::to_string( type.dim() ) + ")";
    }
  }



  // geometryTypeFromString
  // ----------------------

  inline static GeometryType geometryTypeFromString ( const std::string &s )
  {
    typedef GeometryType (*Constructor) ( const std::vector< std::string > & );
    static const char *constructorNames[] = {
        "cube",
        "general",
        "hexahedron",
        "line",
        "none",
        "prism",
        "pyramid",
        "quadrilateral",
        "simplex",
        "tetrahedron",
        "triangle",
        "vertex"
      };
    static const Constructor constructors[]
      = {
          // cube
          [] ( const std::vector< std::string > &args ) {
              if( args.size() != 1 )
                DUNE_THROW( Exception, "GeometryType 'cube' requires integer argument for dimension." );
              return GeometryTypes::cube( std::stoul( args[ 0 ] ) );
            },
          // general
          [] ( const std::vector< std::string > &args ) {
              if( args.size() != 2 )
                DUNE_THROW( Exception, "GeometryType 'general' requires two integer arguments, topologyId and dimension." );
              const auto id = std::stoul( args[ 0 ] );
              const auto dim = std::stoul( args[ 1 ] );
              if( id >= Dune::Impl::numTopologies( dim ) )
                DUNE_THROW( Exception, "Topology id " << id << " too large for dimension " << dim << "." );
              return GeometryType( id, dim );
            },
          // hexahedron
          [] ( const std::vector< std::string > &args ) {
              if( !args.empty() )
                DUNE_THROW( Exception, "GeometryType 'hexahedron' does not require arguments." );
              return GeometryTypes::hexahedron;
            },
          // line
          [] ( const std::vector< std::string > &args ) {
              if( !args.empty() )
                DUNE_THROW( Exception, "GeometryType 'line' does not require arguments." );
              return GeometryTypes::line;
            },
          // none
          [] ( const std::vector< std::string > &args ) {
              if( args.size() != 1 )
                DUNE_THROW( Exception, "GeometryType 'none' requires integer argument for dimension." );
              return GeometryTypes::none( std::stoul( args[ 0 ] ) );
            },
          // prism
          [] ( const std::vector< std::string > &args ) {
              if( !args.empty() )
                DUNE_THROW( Exception, "GeometryType 'prism' does not require arguments." );
              return GeometryTypes::prism;
            },
          // pyramid
          [] ( const std::vector< std::string > &args ) {
              if( !args.empty() )
                DUNE_THROW( Exception, "GeometryType 'pyramid' does not require arguments." );
              return GeometryTypes::pyramid;
            },
          // quadrilateral
          [] ( const std::vector< std::string > &args ) {
              if( !args.empty() )
                DUNE_THROW( Exception, "GeometryType 'quadrilateral' does not require arguments." );
              return GeometryTypes::quadrilateral;
            },
          // simplex
          [] ( const std::vector< std::string > &args ) {
              if( args.size() != 1 )
                DUNE_THROW( Exception, "GeometryType 'simplex' requires integer argument for dimension." );
              return GeometryTypes::simplex( std::stoul( args[ 0 ] ) );
            },
          // tetrahedron
          [] ( const std::vector< std::string > &args ) {
              if( !args.empty() )
                DUNE_THROW( Exception, "GeometryType 'tetrahedron' does not require arguments." );
              return GeometryTypes::tetrahedron;
            },
          // triangle
          [] ( const std::vector< std::string > &args ) {
              if( !args.empty() )
                DUNE_THROW( Exception, "GeometryType 'triangle' does not require arguments." );
              return GeometryTypes::triangle;
            },
          // vertex
          [] ( const std::vector< std::string > &args ) {
              if( !args.empty() )
                DUNE_THROW( Exception, "GeometryType 'vertex' does not require arguments." );
              return GeometryTypes::vertex;
            }
        };
    const std::size_t numConstructors = sizeof( constructorNames ) / sizeof( const char * );

    // find constructor index
    std::size_t n = s.find_first_of( '(' );
    const std::string cName = s.substr( 0, n );
    const std::size_t cIdx = std::lower_bound( constructorNames, constructorNames + numConstructors, cName ) - constructorNames;
    if( (cIdx == numConstructors) || (constructorNames[ cIdx ] != cName) )
      DUNE_THROW( Exception, "No DUNE geometry type constructor named '" << cName << "'." );

    // obtain argument vector
    std::vector< std::string > args;
    if( n != std::string::npos )
    {
      while( s[ n ] != ')' )
      {
        // skip leading spaces
        const std::size_t m = s.find_first_not_of( ' ', n+1 );
        if( m == std::string::npos )
          DUNE_THROW( Exception, "Invalid argument list." );

        // find end of argument
        n = s.find_first_of( ",)", m );
        if( (n == std::string::npos) || (n == m) )
          DUNE_THROW( Exception, "Invalid argument list." );

        // remove trailing spaces from argument
        const std::size_t k = s.find_last_not_of( ' ', n-1 );
        assert( (k != std::string::npos) && (k >= m) );

        args.push_back( s.substr( m, k-m+1 ) );
      }
    }

    // call constructor
    return constructors[ cIdx ]( args );
  }



  namespace Python
  {

    pybind11::class_< GeometryType > registerGeometryType ( pybind11::module scope )
    {
      pybind11::class_< GeometryType > cls( scope, "GeometryType" );

      cls.def( pybind11::init( [] ( const std::string &s ) { return new GeometryType( geometryTypeFromString( s ) ); } ) );

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
      cls.def( "__hash__", [] ( const GeometryType &self ) { return GlobalGeometryTypeIndex::index( self ); } );

      cls.def( "__str__", [] ( const GeometryType &self ) { return to_string( self ); } );

      cls.def_property_readonly( "dim", [] ( const GeometryType &self ) { return self.dim(); } );
      cls.def_property_readonly( "dimension", [] ( const GeometryType &self ) { return self.dim(); } );

      return cls;
    }

  } // namespace Python

} // namespace Dune

#endif // ifndef DUNE_PYTHON_GEOMETRY_TYPE_HH
