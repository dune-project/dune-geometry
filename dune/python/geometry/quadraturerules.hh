// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_PYTHON_GEOMETRY_QUADRATURERULES_HH
#define DUNE_PYTHON_GEOMETRY_QUADRATURERULES_HH

#include <array>
#include <tuple>

#include <dune/common/visibility.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/type.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/numpy.h>


namespace Dune
{

  namespace Python
  {
    template <class Rule>
    auto quadratureToNumpy(const Rule &rule)
    {
      pybind11::array_t< double > p(
          { static_cast< ssize_t >( Rule::d ), static_cast< ssize_t >( rule.size() ) },
          {
            static_cast< ssize_t >( reinterpret_cast< const char * >( &rule[ 0 ].position()[ 1 ] ) - reinterpret_cast< const char * >( &rule[ 0 ].position()[ 0 ] ) ),
            static_cast< ssize_t >( reinterpret_cast< const char * >( &rule[ 1 ].position()[ 0 ] ) - reinterpret_cast< const char * >( &rule[ 0 ].position()[ 0 ] ) )
          },
          &rule[ 0 ].position()[ 0 ]
        );
      pybind11::array_t< double > w(
          { static_cast< ssize_t >( rule.size() ) },
          { static_cast< std::size_t >( reinterpret_cast< const char * >( &rule[ 1 ].weight() ) - reinterpret_cast< const char * >( &rule[ 0 ].weight() ) ) },
          &rule[ 0 ].weight()
        );

      return std::make_pair( p, w );
    }


    template <class Rule>
    auto quadratureToNumpy(pybind11::object self)
    {
      const Rule &rule = pybind11::cast< const Rule & >( self );
      pybind11::array_t< double > p(
          { static_cast< ssize_t >( Rule::d ), static_cast< ssize_t >( rule.size() ) },
          {
            static_cast< ssize_t >( reinterpret_cast< const char * >( &rule[ 0 ].position()[ 1 ] ) - reinterpret_cast< const char * >( &rule[ 0 ].position()[ 0 ] ) ),
            static_cast< ssize_t >( reinterpret_cast< const char * >( &rule[ 1 ].position()[ 0 ] ) - reinterpret_cast< const char * >( &rule[ 0 ].position()[ 0 ] ) )
          },
          &rule[ 0 ].position()[ 0 ],
          self
        );

      pybind11::array_t< double > w(
          { static_cast< ssize_t >( rule.size() ) },
          { static_cast< std::size_t >( reinterpret_cast< const char * >( &rule[ 1 ].weight() ) - reinterpret_cast< const char * >( &rule[ 0 ].weight() ) ) },
          &rule[ 0 ].weight(),
          self
        );

      return std::make_pair( p, w );
    }
    namespace detail
    {

      // registerQuadraturePoint
      // -----------------------

      template< class QP >
      inline void registerQuadraturePoint ( pybind11::object scope, pybind11::class_<QP> cls )
      {
        using pybind11::operator""_a;

        typedef typename QP::Vector Vector;
        typedef typename QP::Field Field;

        cls.def( pybind11::init( [] ( const Vector &x, Field w ) { return new QP( x, w ); } ), "position"_a, "weight"_a );

        cls.def_property_readonly( "position", []( const QP &qp ) -> Vector { return qp.position(); } );
        cls.def_property_readonly( "weight", &QP::weight );

      }



      // registerQuadratureRule
      // ----------------------

      template< class Rule, class... options >
      inline void registerQuadratureRule ( pybind11::object scope,
          pybind11::class_<Rule,options...> cls )
      {
        cls.def_property_readonly( "order", &Rule::order );
        cls.def_property_readonly( "type",  &Rule::type );

        cls.def( "get", [] ( pybind11::object self ) {
            return quadratureToNumpy<Rule>(self);
            } );

        cls.def( "__iter__", [] ( const Rule &rule ) { return pybind11::make_iterator( rule.begin(), rule.end() ); } );
      }

    } // namespace detail



    // registerQuadratureRule
    // ----------------------

    template< class ctype, int dim >
    inline auto registerQuadratureRule ( pybind11::object scope )
    {
      typedef typename Dune::QuadraturePoint< ctype, dim > QP;
      auto quadPointCls = insertClass< QP >( scope, "QuadraturePoint",
            GenerateTypeName("Dune::QuadratePoint",MetaType<ctype>(),dim),
            IncludeFiles{"dune/python/geometry/quadraturerules.hh"});
      if (quadPointCls.second)
        detail::registerQuadraturePoint( scope, quadPointCls.first );

      typedef typename Dune::QuadratureRule< ctype, dim > Rule;
      auto quadRule = insertClass< Rule >(scope, "QuadratureRule" + std::to_string(dim),
            GenerateTypeName("Dune::QuadrateRule",MetaType<ctype>(),dim),
            IncludeFiles{"dune/python/geometry/quadraturerules.hh"});
      if (quadRule.second)
        detail::registerQuadratureRule( scope, quadRule.first );
      return quadRule.first;
    }

    template< class ctype, int ... dim >
    inline auto registerQuadratureRule ( pybind11::object scope, std::integer_sequence< int, dim ... > )
    {
      return std::make_tuple( registerQuadratureRule< ctype >( scope, std::integral_constant< int, dim >() )... );
    }

  } // namespace Python

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_GEOMETRY_QUADRATURERULES_HH
