// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GEOMETRY_QUADRATURERULES_TENSORPRODUCTQUADRATURE_HH
#define DUNE_GEOMETRY_QUADRATURERULES_TENSORPRODUCTQUADRATURE_HH

#ifndef DUNE_INCLUDING_IMPLEMENTATION
#error This is a private header that should not be included directly.
#error Use #include <dune/geometry/quadraturerules.hh> instead.
#endif

#include <algorithm>
#include <bitset>

#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules/jacobiNquadrature.hh>

namespace Dune
{

  /**
   * \brief Quadrature rules constructed by tensor and conical multiplication
   *
   * \tparam ctype Number type used for quadrature point coordinates and weights
   * \tparam dim Dimension of the domain of integration
   */
  template< class ctype, int dim >
  class TensorProductQuadratureRule
    : public QuadratureRule< ctype, dim >
  {
    typedef QuadratureRule<ctype, dim> Base;
    typedef QuadraturePoint<ctype, dim> QPoint;
    typedef typename QPoint::Vector Vector;
    typedef QuadratureRule<ctype,dim-1> BaseQuadrature;

    friend class QuadratureRuleFactory<ctype,dim>;

    TensorProductQuadratureRule (unsigned int topologyId, unsigned int order, QuadratureType::Enum qt)
      : Base( GeometryType(topologyId, dim), order )
    {
      constexpr static int bitSize = sizeof(unsigned int)*8;
      std::bitset<bitSize> baseId(topologyId);
      bool isPrism = baseId[dim-1];
      baseId.reset(dim-1);
      GeometryType baseType(baseId.to_ulong(), dim-1);
      const BaseQuadrature & baseQuad = QuadratureRules<ctype,dim-1>::rule(baseType, order, qt);
      if (isPrism)
        tensorProduct(baseQuad, order, qt);
      else
        conicalProduct(baseQuad, order, qt);
    }

    /**
     * \brief Creates quadrature rule by tensor multiplication of an arbitrary rule with a rule for a one-dimensional domain
     *
     * \param baseQuad Quadrature rule for the base domain
     * \param order Requested order of the one-dimensional rule
     * \param qt Type of the one-dimensional rule
     */
    void tensorProduct(const BaseQuadrature & baseQuad, unsigned int order, QuadratureType::Enum qt)
    {
      typedef QuadratureRule<ctype,1> OneDQuadrature;
      const OneDQuadrature & onedQuad =
        QuadratureRules<ctype,1>::rule(GeometryTypes::line, order, qt);

      const unsigned int baseQuadSize = baseQuad.size();
      for( unsigned int bqi = 0; bqi < baseQuadSize; ++bqi )
      {
        const typename QuadraturePoint<ctype, dim-1>::Vector &
        basePoint = baseQuad[bqi].position( );
        const ctype &baseWeight = baseQuad[bqi].weight( );

        Vector point;
        for( unsigned int i = 0; i < dim-1; ++i )
          point[ i ] = basePoint[ i ];

        const unsigned int onedQuadSize = onedQuad.size();
        for( unsigned int oqi = 0; oqi < onedQuadSize; ++oqi )
        {
          point[ dim-1 ] = onedQuad[oqi].position()[ 0 ];
          this->push_back( QPoint(point, baseWeight * onedQuad[oqi].weight()) );
        }
      }
    }

    /** \brief Creates quadrature rule by conical multiplication of an arbitrary rule with a rule for a one-dimensional domain
     *
     *  This quadrature for \f$B^\circ\f$ is generated from a quadrature for
     *  \f$B\f$ and a 1D quadrature by the so-called Duffy-Transformation
     *  \f$y(x,z) = ((1-z)x,z)^T\f$. Hence, we have
     *  \f[
     *  \int_{B^\circ} f( y )\,\mathrm{d}y
     *  = \int_0^1 \int_B f( (1-z)x, z )\,\mathrm{d}x\,(1-z)^{\dim B}\,\mathrm{d}z.
     *  \f]
     *  Therefore, the 1D quadrature must be at least \f$\dim B\f$ orders higher
     *  than the quadrature for \f$B\f$.
     *
     *  Question: If the polynomials are created via Duffy Transformation, do we
     *            really need a higher quadrature order?
     *
     *  Answer (from OS): Not really.  The official way is to use a Gauss-Jacobi
     *            rule with \f$ \alpha = \dim B, \beta = 0 \f$ for the 1D rule.
     *            That takes care of the term \f$ (1-z)^{\dim B} \f$ without needing
     *            additional orders.  See for example A.H. Stroud, "Approximate Calculation
     *            of Multiple Integrals", Chapters 2.4 and 2.5 for details.
     *            If you want to use plain Gauss-Legendre you do need the additional orders.
     *
     * \note Details on how this construction works can be found in the book by A.H. Stroud,
     *    "Approximate Calculation of Multiple Integrals", Chapter 2.5.
     *
     * \param baseQuad Quadrature rule for the base domain
     * \param order Requested order of the one-dimensional rule
     * \param qt Type of the one-dimensional rule
     */
    void conicalProduct(const BaseQuadrature & baseQuad, unsigned int order, QuadratureType::Enum qt)
    {
      typedef QuadratureRule<ctype,1> OneDQuadrature;

      OneDQuadrature onedQuad;
      bool respectWeightFunction = false;
      if( qt != QuadratureType::GaussJacobi_n_0)
        onedQuad = QuadratureRules<ctype,1>::rule(GeometryTypes::line, order + dim-1, qt);
      else
      {
        onedQuad = JacobiNQuadratureRule1D<ctype>(order,dim-1);
        respectWeightFunction = true;
      }

      const unsigned int baseQuadSize = baseQuad.size();
      for( unsigned int bqi = 0; bqi < baseQuadSize; ++bqi )
      {
        const typename QuadraturePoint<ctype, dim-1>::Vector &
        basePoint = baseQuad[bqi].position( );
        const ctype &baseWeight = baseQuad[bqi].weight( );

        const unsigned int onedQuadSize = onedQuad.size();
        for( unsigned int oqi = 0; oqi < onedQuadSize; ++oqi )
        {
          Vector point;
          point[ dim-1 ] = onedQuad[oqi].position( )[ 0 ];
          const ctype scale = ctype( 1 ) - point[ dim-1 ];
          for( unsigned int i = 0; i < dim-1; ++i )
            point[ i ] = scale * basePoint[ i ];

          ctype weight = baseWeight * onedQuad[oqi].weight( );
          if (!respectWeightFunction)
          {
            for ( unsigned int p = 0; p<dim-1; ++p)
              weight *= scale;                    // pow( scale, dim-1 );
          }

          this->push_back( QPoint(point, weight) );
        }
      }
    }

    static unsigned maxOrder(unsigned int topologyId, QuadratureType::Enum qt)
    {
      constexpr static int bitSize = sizeof(unsigned int)*8;
      std::bitset<bitSize> baseId(topologyId);
      bool isPrism = baseId[dim-1];
      baseId.reset(dim-1);
      GeometryType baseType(baseId.to_ulong(), dim-1);
      unsigned order = QuadratureRules<ctype,dim-1>::maxOrder(baseType, qt);
      if (isPrism)
        order = std::min
          (order, QuadratureRules<ctype,1>::maxOrder(GeometryTypes::line, qt));
      else if (qt != QuadratureType::GaussJacobi_n_0)
        order = std::min
          (order, QuadratureRules<ctype,1>::maxOrder(GeometryTypes::line, qt)-(dim-1));
      else
        order = std::min
          (order, QuadratureRules<ctype,1>::maxOrder(GeometryTypes::line, qt));
      return order;
    }

  };

} // end namspace Dune

#endif // DUNE_GEOMETRY_QUADRATURERULES_TENSORPRODUCTQUADRATURE_HH
