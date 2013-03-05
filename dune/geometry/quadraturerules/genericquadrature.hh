// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_QUADRATURERULES_GENERICQUADRATURE_HH
#define DUNE_GEOMETRY_QUADRATURERULES_GENERICQUADRATURE_HH

#include <vector>

#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/topologyfactory.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/genericgeometry/topologytypes.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // GenericQuadrature
    // -----------------

    /**
     * @brief extends a 1d quadrature to a generic reference elemenet
     *
     * \tparam ctype coordinate type
     * \tparam Topology the topology of the reference element
     *
     * The 1d quadrature must be a std::vector-like container with a constructor
     * taking an order parameter
     **/
    template< class ctype, class Topology >
    class GenericQuadrature;

    /** \brief Generic Quadrature for Point
    **/
    template< class ctype >
    class GenericQuadrature< ctype, Point >
      : public QuadratureRule< ctype, 0 >
    {
      typedef QuadraturePoint< ctype, 0 > QPoint;

    public:
      typedef Point Topology;

      static const unsigned int dimension = 0;
      typedef typename QPoint::Vector Vector;

      GenericQuadrature (unsigned int order, QuadratureType::Enum qt)
        : QuadratureRule<ctype, 0>( GeometryType(Topology::id, dimension), order )
      {
        this->push_back( QPoint( Vector(ctype( 0 ) ), 1) );
      }
    };


    /** \brief Generic Quadrature for Prisms
    **/
    template< class ctype, class BaseTopology >
    class GenericQuadrature< ctype, Prism< BaseTopology > >
      : public QuadratureRule< ctype, Prism< BaseTopology >::dimension >
    {
      typedef QuadratureRule< ctype, Prism< BaseTopology >::dimension > Base;
      typedef QuadraturePoint< ctype, Prism< BaseTopology >::dimension > QPoint;

    public:
      typedef Prism< BaseTopology > Topology;

      static const unsigned int dimension = QPoint::dimension;

      typedef typename QPoint::Vector Vector;

    private:
      typedef GenericQuadrature< ctype, BaseTopology > BaseQuadrature;

    public:
      GenericQuadrature (unsigned int order, QuadratureType::Enum qt)
        : Base( GeometryType(Topology::id, dimension), order )
      {
        typedef QuadratureRule<ctype,1> OneDQuadrature;
        GeometryType onedtype(GeometryType::cube,1);
        const OneDQuadrature & onedQuad =
          QuadratureRules<ctype,1>::rule(onedtype, order, qt);
        BaseQuadrature baseQuad( order, qt );

        const unsigned int baseQuadSize = baseQuad.size();
        for( unsigned int bqi = 0; bqi < baseQuadSize; ++bqi )
        {
          const typename BaseQuadrature::Vector &basePoint = baseQuad[bqi].position( );
          const ctype &baseWeight = baseQuad[bqi].weight( );

          Vector point;
          for( unsigned int i = 0; i < dimension-1; ++i )
            point[ i ] = basePoint[ i ];

          const unsigned int onedQuadSize = onedQuad.size();
          for( unsigned int oqi = 0; oqi < onedQuadSize; ++oqi )
          {
            point[ dimension-1 ] = onedQuad[oqi].position()[ 0 ];
            this->push_back( QPoint(point, baseWeight * onedQuad[oqi].weight()) );
          }
        }
      }
    };


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
    template< class ctype, class BaseTopology >
    class GenericQuadrature< ctype, Pyramid< BaseTopology > >
      : public QuadratureRule< ctype, Pyramid< BaseTopology >::dimension >
    {
      typedef QuadratureRule< ctype, Pyramid< BaseTopology >::dimension > Base;
      typedef QuadraturePoint< ctype, Prism< BaseTopology >::dimension > QPoint;

    public:
      typedef Pyramid< BaseTopology > Topology;

      static const unsigned int dimension = QPoint::dimension;

      typedef typename QPoint::Vector Vector;

    private:
      typedef GenericQuadrature< ctype, BaseTopology > BaseQuadrature;

    public:
      GenericQuadrature (unsigned int order, QuadratureType::Enum qt)
        : Base( GeometryType(Topology::id, dimension), order )
      {
        typedef QuadratureRule<ctype,1> OneDQuadrature;
        GeometryType onedtype(GeometryType::cube,1);
        const OneDQuadrature & onedQuad =
          QuadratureRules<ctype,1>::rule(onedtype, order + dimension-1, qt);
        BaseQuadrature baseQuad( order, qt );

        const unsigned int baseQuadSize = baseQuad.size();
        for( unsigned int bqi = 0; bqi < baseQuadSize; ++bqi )
        {
          const typename BaseQuadrature::Vector &basePoint = baseQuad[bqi].position( );
          const ctype &baseWeight = baseQuad[bqi].weight( );

          const unsigned int onedQuadSize = onedQuad.size();
          for( unsigned int oqi = 0; oqi < onedQuadSize; ++oqi )
          {
            Vector point;
            point[ dimension-1 ] = onedQuad[oqi].position( )[ 0 ];
            const ctype scale = ctype( 1 ) - point[ dimension-1 ];
            for( unsigned int i = 0; i < dimension-1; ++i )
              point[ i ] = scale * basePoint[ i ];

            ctype weight = baseWeight * onedQuad[oqi].weight( );
            for ( unsigned int p = 0; p<dimension-1; ++p)
              weight *= scale;                    // pow( scale, dimension-1 );
            this->push_back( QPoint(point, weight) );
          }
        }
      }
    };

    /**
     * @brief Factory for the generic quadratures
     *
     * This is a Dune::GenericGeometry::TopologyFactory creating
     * GenericQuadrature from a given 1d quadrature
     *
     * \tparam dim dimension of the reference elements contained in the factory
     * \tparam F field in which weight and point of the quadrature are stored
     *
     * Note: the computation of the quadrature points and weights are
     * carried out in the field type of the 1d quadrature which can differ from F.
     **/
    template< int dim, class F >
    struct GenericQuadratureFactory;

    template< int dim, class F >
    struct GenericQuadratureFactoryTraits
    {
      static const unsigned int dimension = dim;
      struct Key
      {
        unsigned int order;
        QuadratureType::Enum qt;
      };
      typedef const QuadratureRule<F,dim> Object;
      typedef GenericQuadratureFactory<dim,F> Factory;
    };

    template< int dim, class F >
    struct GenericQuadratureFactory :
      public TopologyFactory< GenericQuadratureFactoryTraits<dim,F> >
    {
      static const unsigned int dimension = dim;
      typedef F Field;
      typedef GenericQuadratureFactoryTraits<dim,F> Traits;

      typedef typename Traits::Key Key;
      typedef typename Traits::Object Object;

      template< class Topology >
      static Object* createObject ( const Key &key )
      {
        return new Object( GenericQuadrature< F, Topology >( key.order, key.qt ) );
      }
    };


    // GenericQuadratureProvider
    // ---------------------

    /**
     * @brief Singleton factory for the generic quadratures
     *
     * Wrapper for the Dune::GenericGeometry::GenericQuadratureFactory providing
     * singleton storage.
     *
     * \tparam dim dimension of the reference elements contained in the factory
     * \tparam F field in which weight and point of the quadrature are stored
     **/
    template< int dim, class F >
    struct GenericQuadratureProvider
      : public TopologySingletonFactory< GenericQuadratureFactory< dim, F > >
    {};

  }

}

#endif // #ifndef DUNE_GEOMETRY_QUADRATURERULES_GENERICQUADRATURE_HH
