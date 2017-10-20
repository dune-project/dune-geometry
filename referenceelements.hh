#ifndef DUNE_PYTHON_GEOMETRY_REFERENCEELEMENTS_HH
#define DUNE_PYTHON_GEOMETRY_REFERENCEELEMENTS_HH

#include <array>
#include <functional>

#include <dune/common/visibility.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/python/geometry/quadraturerules.hh>
#include <dune/python/pybind11/pybind11.h>

namespace Dune
{

  namespace CorePy
  {

    template< class RefElement, class... options >
    void registerReferenceElements ( pybind11::module module, pybind11::class_< RefElement, options... > cls )
    {
      cls.def( "size", [] ( const RefElement &e, int c ) { return e.size( c ); } );
      cls.def( "size", [] ( const RefElement &e, int i, int c, int cc ) { return e.size( i, c, cc ); } );
      cls.def( "type", [] ( const RefElement &e, int i, int c ) { return e.type( i, c ); } );

      cls.def( "subEntity", []( const RefElement &e, int a,int b,int c,int d) { return e.subEntity(a,b,c,d); } );
      cls.def( "position", &RefElement::position );
      cls.def_property_readonly( "center", [] ( const RefElement &e ) { return e.position( 0, 0 ); } );
      cls.def( "volume", &RefElement::volume );
      cls.def( "integrationOuterNormal", &RefElement::integrationOuterNormal );

      cls.def_property_readonly( "type", [] ( const RefElement &e ) { return e.type(); } );

      registerQuadratureRule<typename RefElement::ctype, RefElement::dimension>( module );

      module.def( "general", [] ( const GeometryType &gt ) {
            return Dune::ReferenceElements< typename RefElement::ctype, RefElement::dimension >::general( gt );
      }, pybind11::return_value_policy::reference );
      module.def( "rule", [] (const GeometryType &gt, int order) {
            return Dune::QuadratureRules< typename RefElement::ctype, RefElement::dimension >::rule( gt, order );
      }, pybind11::return_value_policy::reference );
    }

  } // namespace CorePy

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_GEOMETRY_REFERENCEELEMENTS_HH
