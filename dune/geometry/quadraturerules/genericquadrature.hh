// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_QUADRATURERULES_GENERICQUADRATURE_HH
#define DUNE_GEOMETRY_QUADRATURERULES_GENERICQUADRATURE_HH

#include <bitset>

namespace Dune
{

  /**
   * \brief Generic Quadrature rules using prism and pyramid construction (dynamic implementation)
   */
  template< class ctype, int dim >
  class GenericQuadratureRule
    : public QuadratureRule< ctype, dim >
  {
    typedef QuadratureRule<ctype, dim> Base;
    typedef QuadraturePoint<ctype, dim> QPoint;
    typedef typename QPoint::Vector Vector;
    typedef QuadratureRule<ctype,dim-1> BaseQuadrature;

    friend class QuadratureRuleFactory<ctype,dim>;

    GenericQuadratureRule (unsigned int topologyId, unsigned int order, QuadratureType::Enum qt)
      : Base( GeometryType(topologyId, dim), order )
    {
      enum { bitSize = sizeof(unsigned int)*8 };
      std::bitset<bitSize> baseId(topologyId);
      bool isPrism = baseId[dim-1];
      baseId.reset(dim-1);
      GeometryType baseType(baseId.to_ulong(), dim-1);
      const BaseQuadrature & baseQuad = QuadratureRules<ctype,dim-1>::rule(baseType, order, qt);
      if (isPrism)
        create_prism(baseQuad, order, qt);
      else
        create_pyramid(baseQuad, order, qt);
    }

    /**
     * \brief Generic Quadrature for Prisms
     */
    void create_prism(const BaseQuadrature & baseQuad, unsigned int order, QuadratureType::Enum qt)
    {
      typedef QuadratureRule<ctype,1> OneDQuadrature;
      GeometryType onedType(GeometryType::cube,1);
      const OneDQuadrature & onedQuad =
        QuadratureRules<ctype,1>::rule(onedType, order, qt);

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

    /** \brief Generic Quadrature for Pyramids
     *
     *  This quadrature for \f$B^\circ\f$ is generated from a quadrature for
     *  \f$B\f$ and a 1D quadrature by the so-called Duffy-Transformation
     *  \f$y(x,z) = ((1-z)x,y)^T\f$. Hence, we have
     *  \f[
     *  \int_{B^\circ} f( y )\,\mathrm{d}y
     *  = \int_0^1 \int_B f( (1-z)x, z )\,\mathrm{d}x\,(1-z)^{\dim B}\,\mathrm{d}z.
     *  \f]
     *  Therefore, the 1D quadrature must be at least \f$\dim B\f$ orders higher
     *  than the quadrature for \f$B\f$.
     *
     *  Question: If the polynomials are created via Duffy Transformation, do we
     *            really need a higher quadrature order?
     */
    void create_pyramid(const BaseQuadrature & baseQuad, unsigned int order, QuadratureType::Enum qt)
    {
      typedef QuadratureRule<ctype,1> OneDQuadrature;
      GeometryType onedtype(GeometryType::cube,1);
      const OneDQuadrature & onedQuad =
        QuadratureRules<ctype,1>::rule(onedtype, order + dim-1, qt);

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
          for ( unsigned int p = 0; p<dim-1; ++p)
            weight *= scale;                    // pow( scale, dim-1 );
          this->push_back( QPoint(point, weight) );
        }
      }
    }
  };

}

#endif // #ifndef DUNE_GEOMETRY_QUADRATURERULES_GENERICQUADRATURE_HH
