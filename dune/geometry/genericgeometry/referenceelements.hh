// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GEOMETRY_GENERICGEOMETRY_REFERENCEELEMENTS_HH
#define DUNE_GEOMETRY_GENERICGEOMETRY_REFERENCEELEMENTS_HH

/** \file
 *  \brief Implements some reference element functionality needed by the generic geometries
 *  \warning This is an internal header.  Do not include it from outside of dune-geometry.
 */
#include <dune/common/array.hh>
#include <dune/common/fvector.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/visibility.hh>

#include <dune/geometry/genericgeometry/referencedomain.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // ReferenceElement
    // ----------------

    template< class Topology, class ctype >
    struct ReferenceElement
    {
      static const unsigned int topologyId = Topology :: id;
      static const unsigned int dimension = Topology :: dimension;

      static const unsigned int numCorners = Topology :: numCorners;
      static const unsigned int numNormals = ReferenceDomain< Topology > :: numNormals;

      typedef FieldVector< ctype, dimension > CoordinateType;

      template< unsigned int codim >
      struct Codim
      {
        enum { size = Size< Topology, codim > :: value };
      };

      template< unsigned int codim, unsigned int subcodim >
      static unsigned int subNumbering ( unsigned int i, unsigned int j )
      {
        return SubTopologyNumbering< Topology, codim, subcodim > :: number( i, j );
      }

      template< unsigned int codim, unsigned int subcodim >
      static unsigned int size ( unsigned int i )
      {
        return SubTopologySize< Topology, codim, subcodim > :: size( i );
      }

      /** \brief Return the element barycenter
       */
      static const FieldVector< ctype, dimension > &
      baryCenter ()
      {
        return instance().baryCenter_;
      }

      static const CoordinateType &corner ( unsigned int i )
      {
        assert( i < numCorners );
        return instance().corners_[ i ];
      }

      static bool checkInside ( const CoordinateType &x )
      {
        return ReferenceDomain< Topology >::checkInside( x );
      }

      static const CoordinateType &
      integrationOuterNormal ( unsigned int i )
      {
        assert( i < numNormals );
        return instance().normals_[ i ];
      }

      static ctype volume ()
      {
        return ReferenceDomain< Topology > :: template volume< ctype >();
      }

      DUNE_EXPORT static const ReferenceElement &instance ()
      {
        static ReferenceElement inst;
        return inst;
      }

    private:
      class BaryCenterArray;

      ReferenceElement ()
      {
        for( unsigned int i = 0; i < numCorners; ++i )
          ReferenceDomain< Topology > :: corner( i, corners_[ i ] );
        for( unsigned int i = 0; i < numNormals; ++i )
          ReferenceDomain< Topology > :: integrationOuterNormal( i, normals_[ i ] );

        // Compute the element barycenter
        typedef SubTopologyNumbering< Topology, 0, dimension > Numbering;
        typedef SubTopologySize< Topology, 0, dimension > Size;

        baryCenter_ = 0;
        const unsigned int numCorners = Size :: size( 0 );
        for( unsigned int k = 0; k < numCorners; ++k )
        {
          unsigned int j = Numbering :: number( 0, k );

          CoordinateType y;
          ReferenceDomain< Topology > :: corner( j, y );
          baryCenter_ += y;
        }
        baryCenter_ *= ctype( 1 ) / ctype( numCorners );
      }

      Dune::array< CoordinateType, numCorners > corners_;
      CoordinateType baryCenter_;
      Dune::array< CoordinateType, numNormals > normals_;
    };

  }

}

#endif // DUNE_GEOMETRY_GENERICGEOMETRY_REFERENCEELEMENTS_HH
