// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_PYTHON_GEOMETRY_REFERENCEELEMENTS_HH
#define DUNE_PYTHON_GEOMETRY_REFERENCEELEMENTS_HH

#include <array>
#include <functional>
#include <string>

#include <dune/common/visibility.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/python/common/fvecmatregistry.hh>

#include <dune/python/geometry/quadraturerules.hh>
#include <dune/python/pybind11/pybind11.h>

namespace Dune
{

  namespace Python
  {

    namespace detail
    {

      // referenceElementSize
      // --------------------

      template< class RefElement >
      inline static int referenceElementSize ( const RefElement &refElement, int c )
      {
        if( (c < 0) || (c > RefElement::dimension) )
          throw pybind11::value_error( "Invalid codimension: " + std::to_string( c ) + " (must be in [0, " + std::to_string( RefElement::dimension ) + "])." );
        return refElement.size( c );
      }

      template< class RefElement >
      inline static int referenceElementSize ( const RefElement &refElement, int i, int c, int cc )
      {
        const int size = detail::referenceElementSize( refElement, c );
        if( (i < 0) || (i >= size) )
          throw pybind11::value_error( "Invalid index: " + std::to_string( i ) + " (must be in [0, " + std::to_string( size ) + "))." );
        if( (cc < c) || (cc > RefElement::dimension) )
          throw pybind11::value_error( "Invalid codimension: " + std::to_string( cc ) + " (must be in [" + std::to_string( c ) + ", " + std::to_string( RefElement::dimension ) + "])." );
        return refElement.size( i, c, cc );
      }



      // referenceElementSubEntity
      // -------------------------

      template< class RefElement >
      inline static int referenceElementSubEntity ( const RefElement &refElement, int i, int c, int ii, int cc )
      {
        const int size = detail::referenceElementSize( refElement, i, c, cc );
        if( (ii < 0) || (ii >= size) )
          throw pybind11::value_error( "Invalid index: " + std::to_string( i ) + " (must be in [0, " + std::to_string( size ) + "))." );
        return refElement.subEntity( i, c, ii, cc );
      }



      // referenceElementPosition
      // ------------------------

      template< class RefElement >
      inline static auto referenceElementPosition ( const RefElement &refElement, int i, int c )
      {
        const int size = detail::referenceElementSize( refElement, c );
        if( (i < 0) || (i >= size) )
          throw pybind11::value_error( "Invalid index: " + std::to_string( i ) + " (must be in [0, " + std::to_string( size ) + "))." );
        return refElement.position( i, c );
      }


      // referenceElementType
      // --------------------

      template< class RefElement >
      inline static GeometryType referenceElementType ( const RefElement &refElement, int i, int c )
      {
        const int size = detail::referenceElementSize( refElement, c );
        if( (i < 0) || (i >= size) )
          throw pybind11::value_error( "Invalid index: " + std::to_string( i ) + " (must be in [0, " + std::to_string( size ) + "))." );
        return refElement.type( i, c );
      }

    } // namespace detail


    template< class RefElement, class... options >
    void registerReferenceElements ( pybind11::module module, pybind11::class_< RefElement, options... > cls )
    {
      using pybind11::operator""_a;

      pybind11::options opts;
      opts.disable_function_signatures();

      static const std::size_t dimension = RefElement::dimension;
      typedef typename RefElement::ctype ctype;
      registerFieldVecMat<FieldVector<ctype,dimension>>::apply();
      registerFieldVecMat<FieldMatrix<ctype,dimension,dimension>>::apply();

      cls.def_property_readonly( "dimension", [] ( pybind11::object ) { return pybind11::int_( RefElement::dimension ); } );

      cls.def( "size", [] ( const RefElement &self, int c ) {
          return detail::referenceElementSize( self, c );
        }, "codim"_a );
      cls.def( "size", [] ( const RefElement &self, int i, int c, int cc ) {
          return detail::referenceElementSize( self, i, c, cc );
         } );
      cls.def( "size", [] ( const RefElement &self, std::tuple< int, int > e, int cc ) {
          return detail::referenceElementSize( self, std::get< 0 >( e ), std::get< 1 >( e ), cc );
        }, "subentity"_a, "codim"_a );

      cls.def( "type", [] ( const RefElement &self, int i, int c ) {
          return detail::referenceElementType( self, i, c );
        }, "index"_a, "codim"_a );
      cls.def( "type", [] ( const RefElement &self, std::tuple< int, int > e ) {
          return detail::referenceElementType( self, std::get< 0 >( e ), std::get< 1 >( e ) );
        }, "subentity"_a );
      cls.def( "types", [] ( const RefElement &self, int c ) {
          const int size = detail::referenceElementSize( self, c );
          pybind11::tuple types( size );
          for( int i = 0; i < size; ++i )
            types[ i ] = pybind11::cast( self.type( i, c ) );
          return types;
        }, "codim"_a,
        R"doc(
          get geometry types of subentities

          Args:
              codim:    codimension of subentities

          Returns:
              tuple containing geometry type for each subentity
        )doc" );

      cls.def( "subEntity", [] ( const RefElement &self, int i, int c, int ii, int cc ) {
          return detail::referenceElementSubEntity( self, i, c, ii, cc );
        } );
      cls.def( "subEntity", [] ( const RefElement &self, std::tuple< int, int > e, std::tuple< int, int > ee ) {
          return detail::referenceElementSubEntity( self, std::get< 0 >( e ), std::get< 1 >( e ), std::get< 0 >( ee ), std::get< 1 >( ee ) );
        } );
      cls.def( "subEntities", [] ( const RefElement &self, int i, int c, int cc ) {
          const int size = detail::referenceElementSize( self, i, c, cc );
          pybind11::tuple subEntities( size );
          for( int ii = 0; ii < size; ++ii )
            subEntities[ ii ] = pybind11::int_( self.subEntity( i, c, ii, cc ) );
          return subEntities;
        } );
      cls.def( "subEntities", [] ( const RefElement &self, std::tuple< int, int > e, int cc ) {
          const int size = detail::referenceElementSize( self, std::get< 0 >( e ), std::get< 1 >( e ), cc );
          pybind11::tuple subEntities( size );
          for( int ii = 0; ii < size; ++ii )
            subEntities[ ii ] = pybind11::int_( self.subEntity( std::get< 0 >( e ), std::get< 1 >( e ), ii, cc ) );
          return subEntities;
        }, "subentity"_a, "codim"_a );

      cls.def( "position", [] ( const RefElement &self, int i, int c ) {
          return detail::referenceElementPosition( self, i, c );
        }, "index"_a, "codim"_a );
      cls.def( "position", [] ( const RefElement &self, std::tuple< int, int > e ) {
          return detail::referenceElementPosition( self, std::get< 0 >( e ), std::get< 1 >( e ) );
        }, "subentity"_a );
      cls.def( "positions", [] ( const RefElement &self, int c ) {
          const int size = detail::referenceElementSize( self, c );
          pybind11::tuple positions( size );
          for( int i = 0; i < size; ++i )
            positions[ i ] = pybind11::cast( self.position( i, c ) );
          return positions;
        }, "codim"_a,
        R"doc(
          get barycenters of subentities

          Args:
              codim:    codimension of subentities

          Returns:
              tuple containing barycenter for each subentity
        )doc" );

#if 0
      // Bug: This property overwrite the method "type"
      cls.def_property_readonly( "type", [] ( const RefElement &self ) {
          return self.type();
        },
        R"doc(
          geometry type of reference element
        )doc" );
#endif
      cls.def_property_readonly( "center", [] ( const RefElement &self ) {
          return self.position( 0, 0 );
        },
        R"doc(
          barycenter of reference domain
        )doc" );
      cls.def_property_readonly( "corners", [] ( const RefElement &self ) {
          const int size = self.size( RefElement::dimension );
          pybind11::tuple corners( size );
          for( int i = 0; i < size; ++i )
            corners[ i ] = pybind11::cast( self.position( i, RefElement::dimension ) );
          return corners;
        },
        R"doc(
          corners of reference domain
        )doc" );
      cls.def_property_readonly( "volume", [] ( const RefElement &self ) {
          return self.volume();
        },
        R"doc(
          volume of reference domain
        )doc" );

      if( RefElement::dimension > 0 )
      {
        cls.def( "integrationOuterNormal", [] ( const RefElement &self, int i ) {
            if( (i < 0) || (i >= self.size( 1 )) )
              throw pybind11::value_error( "Invalid index: " + std::to_string( i ) + " (must be in [0, " + std::to_string( self.size( 1 ) ) + "))." );
            return self.integrationOuterNormal( i );
          }, "index"_a );
        cls.def_property_readonly( "integrationOuterNormals", [] ( const RefElement &self ) {
            const int size = self.size( 1 );
            pybind11::tuple integrationOuterNormals( size );
            for( int i = 0; i < size; ++i )
              integrationOuterNormals[ i ] = pybind11::cast( self.integrationOuterNormal( i ) );
            return integrationOuterNormals;
          } );
      }

      registerQuadratureRule<typename RefElement::ctype, RefElement::dimension>( module );

      module.def( "general", [] ( const GeometryType &gt ) -> pybind11::object {
          if( gt.isNone() )
            return pybind11::none();
          else
            return pybind11::cast( Dune::ReferenceElements< typename RefElement::ctype, RefElement::dimension >::general( gt ), pybind11::return_value_policy::reference );
        } );
      module.def( "rule", [] (const GeometryType &gt, int order) {
            return Dune::QuadratureRules< typename RefElement::ctype, RefElement::dimension >::rule( gt, order );
      }, pybind11::return_value_policy::reference );
    }

  } // namespace Python

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_GEOMETRY_REFERENCEELEMENTS_HH
