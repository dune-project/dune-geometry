// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_QUADRATURERULES_GAUSSQUADRATURE_HH
#define DUNE_GEOMETRY_QUADRATURERULES_GAUSSQUADRATURE_HH

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/quadraturerules/genericquadrature.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // GaussPoints
    // -----------

    /**
     * @brief Gauss quadrature points and weights in 1d.
     **/
    template< class F>
    class GaussPoints
      : public std::vector< QuadraturePoint<F,1> >
    {
      typedef std::vector< QuadraturePoint<F,1> > Base;
      typedef GaussPoints< F > This;

    public:
      typedef F Field;
      // n is number of points required
      explicit GaussPoints ( unsigned int n )
      {
        Base::reserve( n );
        const QuadratureRule<Field,1>& points =
          QuadratureRules<Field,1>::rule(GeometryType(GeometryType::cube,1),
                                         2*n-2, QuadratureType::Gauss);
        for( unsigned int i = 0; i < n; ++i )
        {
          QuadraturePoint<Field,1> q( points[i].position()[0], points[i].weight() );
          Base::push_back( q );
        }
      }
      // order is the maximal order of polynomials which can be exactly integrated
      static unsigned int minPoints( unsigned int order )
      {
        return (order+2)/2;
      }
    };

    /**
     * @brief Singleton provider for Gauss quadratures
     *
     * \tparam dim dimension of the reference elements contained in the factory
     * \tparam F field in which weight and point of the quadrature are stored
     * \tparam CF the compute field for the points and weights
     **/
    template< int dim, class F, class CF=F >
    struct GaussQuadratureProvider
      : public TopologySingletonFactory< GenericQuadratureFactory< dim, F, GaussPoints<CF> > >
    {};

  }

}

#endif // #ifndef DUNE_GEOMETRY_QUADRATURERULES_GAUSSQUADRATURE_HH
