// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_GENERICGEOMETRY_REFERENCEDOMAIN_HH
#define DUNE_GEOMETRY_GENERICGEOMETRY_REFERENCEDOMAIN_HH

#include <algorithm>

#include <dune/common/array.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/genericgeometry/topologytypes.hh>
#include <dune/geometry/genericgeometry/subtopologies.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class Topology >
    struct ReferenceDomain;



    // ReferenceDomain
    // ---------------

    template< class Topology >
    class ReferenceDomainBase;

    /** \cond */
    template<>
    class ReferenceDomainBase< Point >
    {
      typedef Point Topology;

      friend struct ReferenceDomain< Topology >;
      friend class ReferenceDomainBase< Prism< Topology > >;
      friend class ReferenceDomainBase< Pyramid< Topology > >;

      static const unsigned int numNormals = 0;

      template< class ctype, int dim >
      static void corner ( unsigned int i, FieldVector< ctype, dim > &n )
      {
        assert( i < Topology::numCorners );
      }

      template< class ctype, int dim >
      static bool
      checkInside ( const FieldVector< ctype, dim > &x, ctype factor )
      {
        return true;
      }

      template< class ctype, int dim >
      static void
      integrationOuterNormal ( unsigned int i, FieldVector< ctype, dim > &n )
      {
        assert( i < numNormals );
      }

      template< class ctype >
      static ctype volume ()
      {
        return ctype( 1 );
      }
    };


    template< class BaseTopology >
    class ReferenceDomainBase< Prism< BaseTopology > >
    {
      typedef Prism< BaseTopology > Topology;

      friend struct ReferenceDomain< Topology >;
      friend class ReferenceDomainBase< Prism< Topology > >;
      friend class ReferenceDomainBase< Pyramid< Topology > >;

      static const unsigned int numNormals = Size< Topology, 1 >::value;

      static const unsigned int dimension = Topology::dimension;
      static const unsigned int myindex = dimension - 1;

      template< class ctype, int dim >
      static void corner ( unsigned int i, FieldVector< ctype, dim > &x )
      {
        assert( i < Topology::numCorners );
        const unsigned int j = i % BaseTopology::numCorners;
        ReferenceDomainBase< BaseTopology >::corner( j, x );
        if( i >= BaseTopology::numCorners )
          x[ myindex ] = ctype( 1 );
      }

      template< class ctype, int dim >
      static bool
      checkInside ( const FieldVector< ctype, dim > &x, ctype factor )
      {
        const ctype xn = x[ myindex ];
        const ctype cxn = factor - xn;
        return (xn > -1e-12) && (cxn > -1e-12)
               && ReferenceDomainBase< BaseTopology >::checkInside( x, factor );
      }

      template< class ctype, int dim >
      static void
      integrationOuterNormal ( unsigned int i, FieldVector< ctype, dim > &n )
      {
        typedef ReferenceDomainBase< BaseTopology > BaseReferenceDomain;

        if( i >= BaseReferenceDomain::numNormals )
        {
          const unsigned int j = i - BaseReferenceDomain::numNormals;
          n[ myindex ] = (j == 0 ? ctype( -1 ) : ctype( 1 ));
        }
        else
          BaseReferenceDomain::integrationOuterNormal( i, n );
      }

      template< class ctype >
      static ctype volume ()
      {
        typedef ReferenceDomainBase< BaseTopology > BaseReferenceDomain;
        return BaseReferenceDomain::template volume< ctype >();
      }
    };


    template< class BaseTopology >
    class ReferenceDomainBase< Pyramid< BaseTopology > >
    {
      typedef Pyramid< BaseTopology > Topology;

      friend struct ReferenceDomain< Topology >;
      friend class ReferenceDomainBase< Prism< Topology > >;
      friend class ReferenceDomainBase< Pyramid< Topology > >;

      static const unsigned int numNormals = Size< Topology, 1 >::value;

      static const unsigned int dimension = Topology::dimension;
      static const unsigned int myindex = dimension - 1;

      template< bool >
      struct MultiDimensional
      {
        template< class ctype, int dim >
        static void
        integrationOuterNormal ( unsigned int i, FieldVector< ctype, dim > &n )
        {
          multiDimensionalIntegrationOuterNormal( i, n );
        }
      };

      template< bool >
      struct OneDimensional
      {
        template< class ctype, int dim >
        static void
        integrationOuterNormal ( unsigned int i, FieldVector< ctype, dim > &n )
        {
          n[ myindex ] = (i > 0) ? ctype( 1 ) : ctype( -1 );
        }
      };

      template< class ctype, int dim >
      static void corner ( unsigned int i, FieldVector< ctype, dim > &x )
      {
        assert( i < Topology::numCorners );
        if( i < BaseTopology::numCorners )
          ReferenceDomainBase< BaseTopology >::corner( i, x );
        else
          x[ myindex ] = ctype( 1 );
      }

      template< class ctype, int dim >
      static bool
      checkInside ( const FieldVector< ctype, dim > &x, ctype factor )
      {
        const ctype xn = x[ myindex ];
        const ctype cxn = factor - xn;
        return (xn > -1e-12) && (cxn > -1e-12)
               && ReferenceDomainBase< BaseTopology >::checkInside( x, cxn );
      }

      template< class ctype, int dim >
      static void
      multiDimensionalIntegrationOuterNormal ( unsigned int i, FieldVector< ctype, dim > &n )
      {
        typedef ReferenceDomainBase< BaseTopology > BaseReferenceDomain;
        typedef SubTopologyNumbering< BaseTopology, 1, dimension-2 > Numbering;

        if( i > 0 )
        {
          const unsigned int j = Numbering::number( i-1, 0 );
          FieldVector< ctype, dim > x( ctype( 0 ) );
          BaseReferenceDomain::corner( j, x );

          BaseReferenceDomain::integrationOuterNormal ( i-1, n );
          n[ myindex ] = (x * n);
        }
        else
          n[ myindex ] = ctype( -1 );
      }

      template< class ctype, int dim >
      static void
      integrationOuterNormal ( unsigned int i, FieldVector< ctype, dim > &n )
      {
        conditional< (dimension > 1), MultiDimensional<true>, OneDimensional<false> > :: type
        ::integrationOuterNormal( i, n );
      }

      template< class ctype >
      static ctype volume ()
      {
        typedef ReferenceDomainBase< BaseTopology > BaseReferenceDomain;
        const ctype baseVolume = BaseReferenceDomain::template volume< ctype >();
        return baseVolume / ctype( (unsigned int)(dimension) ); // linker problem when using dimension directly
      }
    };
    /** \endcond */



    // ReferenceDomain
    // ---------------

    template< class Topology >
    struct ReferenceDomain
    {
      static const unsigned int numCorners = Topology::numCorners;
      static const unsigned int dimension = Topology::dimension;

      static const unsigned int numNormals
        = ReferenceDomainBase< Topology >::numNormals;

      template< class ctype >
      static void corner ( unsigned int i, FieldVector< ctype, dimension > &x )
      {
        x = ctype( 0 );
        ReferenceDomainBase< Topology >::corner( i, x );
      }

      template< class ctype >
      static bool checkInside ( const FieldVector< ctype, dimension > &x )
      {
        return ReferenceDomainBase< Topology >::checkInside( x, ctype( 1 ) );
      }

      template< class ctype >
      static void
      integrationOuterNormal ( unsigned int i, FieldVector< ctype, dimension > &n )
      {
        n = ctype( 0 );
        return ReferenceDomainBase< Topology >::integrationOuterNormal( i, n );
      }

      template< class ctype >
      static ctype volume ()
      {
        return ReferenceDomainBase< Topology >::template volume< ctype >();
      }
    };



    // checkInside
    // -----------

    template< class ct, int cdim >
    inline bool
    checkInside ( unsigned int topologyId, int dim, const FieldVector< ct, cdim > &x, ct tolerance, ct factor = ct( 1 ) )
    {
      assert( (dim >= 0) && (dim <= cdim) );
      assert( topologyId < numTopologies( dim ) );

      if( dim > 0 )
      {
        const ct baseFactor = (isPrism( topologyId, dim ) ? factor : factor - x[ dim-1 ]);
        if( (x[ dim-1 ] > -tolerance) && (factor - x[ dim-1 ] > -tolerance) )
          return checkInside< ct, cdim >( baseTopologyId( topologyId, dim ), dim-1, x, tolerance, baseFactor );
        else
          return false;
      }
      else
        return true;
    }



    // referenceCorners
    // ----------------

    template< class ct, int cdim >
    inline unsigned int
    referenceCorners ( unsigned int topologyId, int dim, FieldVector< ct, cdim > *corners )
    {
      assert( (dim >= 0) && (dim <= cdim) );
      assert( topologyId < numTopologies( dim ) );

      if( dim > 0 )
      {
        const unsigned int nBaseCorners
          = referenceCorners( baseTopologyId( topologyId, dim ), dim-1, corners );
        assert( nBaseCorners == size( baseTopologyId( topologyId, dim ), dim-1, dim-1 ) );
        if( isPrism( topologyId, dim ) )
        {
          std::copy( corners, corners + nBaseCorners, corners + nBaseCorners );
          for( unsigned int i = 0; i < nBaseCorners; ++i )
            corners[ i+nBaseCorners ][ dim-1 ] = ct( 1 );
          return 2*nBaseCorners;
        }
        else
        {
          corners[ nBaseCorners ] = FieldVector< ct, cdim >( ct( 0 ) );
          corners[ nBaseCorners ][ dim-1 ] = ct( 1 );
          return nBaseCorners+1;
        }
      }
      else
      {
        *corners = FieldVector< ct, cdim >( ct( 0 ) );
        return 1;
      }
    }



    // referenceVolume
    // ---------------

    unsigned long referenceVolumeInverse ( unsigned int topologyId, int dim );

    template< class ct >
    inline ct referenceVolume ( unsigned int topologyId, int dim )
    {
      return ct( 1 ) / ct( referenceVolumeInverse( topologyId, dim ) );
    }



    // referenceOrigins
    // ----------------

    template< class ct, int cdim >
    inline unsigned int
    referenceOrigins ( unsigned int topologyId, int dim, int codim, FieldVector< ct, cdim > *origins )
    {
      assert( (dim >= 0) && (dim <= cdim) );
      assert( topologyId < numTopologies( dim ) );
      assert( (codim >= 0) && (codim <= dim) );

      if( codim > 0 )
      {
        const unsigned int baseId = baseTopologyId( topologyId, dim );
        if( isPrism( topologyId, dim ) )
        {
          const unsigned int n = (codim < dim ? referenceOrigins( baseId, dim-1, codim, origins ) : 0);
          const unsigned int m = referenceOrigins( baseId, dim-1, codim-1, origins+n );
          for( unsigned int i = 0; i < m; ++i )
          {
            origins[ n+m+i ] = origins[ n+i ];
            origins[ n+m+i ][ dim-1 ] = ct( 1 );
          }
          return n+2*m;
        }
        else
        {
          const unsigned int m = referenceOrigins( baseId, dim-1, codim-1, origins );
          if( codim == dim )
          {
            origins[ m ] = FieldVector< ct, cdim >( ct( 0 ) );
            origins[ m ][ dim-1 ] = ct( 1 );
            return m+1;
          }
          else
            return m+referenceOrigins( baseId, dim-1, codim, origins+m );
        }
      }
      else
      {
        origins[ 0 ] = FieldVector< ct, cdim >( ct( 0 ) );
        return 1;
      }
    }



    // referenceEmbeddings
    // -------------------

    template< class ct, int cdim, int mydim >
    inline unsigned int
    referenceEmbeddings ( unsigned int topologyId, int dim, int codim,
                          FieldVector< ct, cdim > *origins,
                          FieldMatrix< ct, mydim, cdim > *jacobianTransposeds )
    {
      assert( (0 <= codim) && (codim <= dim) && (dim <= cdim) );
      assert( (dim - codim <= mydim) && (mydim <= cdim) );
      assert( topologyId < numTopologies( dim ) );

      if( codim > 0 )
      {
        const unsigned int baseId = baseTopologyId( topologyId, dim );
        if( isPrism( topologyId, dim ) )
        {
          const unsigned int n = (codim < dim ? referenceEmbeddings( baseId, dim-1, codim, origins, jacobianTransposeds ) : 0);
          for( unsigned int i = 0; i < n; ++i )
            jacobianTransposeds[ i ][ dim-codim-1 ][ dim-1 ] = ct( 1 );

          const unsigned int m = referenceEmbeddings( baseId, dim-1, codim-1, origins+n, jacobianTransposeds+n );
          std::copy( origins+n, origins+n+m, origins+n+m );
          std::copy( jacobianTransposeds+n, jacobianTransposeds+n+m, jacobianTransposeds+n+m );
          for( unsigned int i = 0; i < m; ++i )
            origins[ n+m+i ][ dim-1 ] = ct( 1 );

          return n+2*m;
        }
        else
        {
          const unsigned int m = referenceEmbeddings( baseId, dim-1, codim-1, origins, jacobianTransposeds );
          if( codim == dim )
          {
            origins[ m ] = FieldVector< ct, cdim >( ct( 0 ) );
            origins[ m ][ dim-1 ] = ct( 1 );
            jacobianTransposeds[ m ] = FieldMatrix< ct, mydim, cdim >( ct( 0 ) );
            return m+1;
          }
          else
          {
            const unsigned int n = referenceEmbeddings( baseId, dim-1, codim, origins+m, jacobianTransposeds+m );
            for( unsigned int i = 0; i < n; ++i )
            {
              for( int k = 0; k < dim-1; ++k )
                jacobianTransposeds[ m+i ][ dim-codim-1 ][ k ] = -origins[ m+i ][ k ];
              jacobianTransposeds[ m+i ][ dim-codim-1 ][ dim-1 ] = ct( 1 );
            }
            return m+n;
          }
        }
      }
      else
      {
        origins[ 0 ] = FieldVector< ct, cdim >( ct( 0 ) );
        jacobianTransposeds[ 0 ] = FieldMatrix< ct, mydim, cdim >( ct( 0 ) );
        for( int k = 0; k < dim; ++k )
          jacobianTransposeds[ 0 ][ k ][ k ] = ct( 1 );
        return 1;
      }
    }



    // referenceIntegrationOuterNormals
    // --------------------------------

    template< class ct, int cdim >
    inline unsigned int
    referenceIntegrationOuterNormals ( unsigned int topologyId, int dim,
                                       const FieldVector< ct, cdim > *origins,
                                       FieldVector< ct, cdim > *normals )
    {
      assert( (dim > 0) && (dim <= cdim) );
      assert( topologyId < numTopologies( dim ) );

      if( dim > 1 )
      {
        const unsigned int baseId = baseTopologyId( topologyId, dim );
        if( isPrism( topologyId, dim ) )
        {
          const unsigned int numBaseFaces
            = referenceIntegrationOuterNormals( baseId, dim-1, origins, normals );

          for( unsigned int i = 0; i < 2; ++i )
          {
            normals[ numBaseFaces+i ] = FieldVector< ct, cdim >( ct( 0 ) );
            normals[ numBaseFaces+i ][ dim-1 ] = ct( 2*int( i )-1 );
          }

          return numBaseFaces+2;
        }
        else
        {
          normals[ 0 ] = FieldVector< ct, cdim >( ct( 0 ) );
          normals[ 0 ][ dim-1 ] = ct( -1 );

          const unsigned int numBaseFaces
            = referenceIntegrationOuterNormals( baseId, dim-1, origins+1, normals+1 );
          for( unsigned int i = 1; i <= numBaseFaces; ++i )
            normals[ i ][ dim-1 ] = normals[ i ]*origins[ i ];

          return numBaseFaces+1;
        }
      }
      else
      {
        for( unsigned int i = 0; i < 2; ++i )
        {
          normals[ i ] = FieldVector< ct, cdim >( ct( 0 ) );
          normals[ i ][ 0 ] = ct( 2*int( i )-1 );
        }

        return 2;
      }
    }

    template< class ct, int cdim >
    inline unsigned int
    referenceIntegrationOuterNormals ( unsigned int topologyId, int dim,
                                       FieldVector< ct, cdim > *normals )
    {
      assert( (dim > 0) && (dim <= cdim) );

      FieldVector< ct, cdim > *origins
        = new FieldVector< ct, cdim >[ size( topologyId, dim, 1 ) ];
      referenceOrigins( topologyId, dim, 1, origins );

      const unsigned int numFaces
        = referenceIntegrationOuterNormals( topologyId, dim, origins, normals );
      assert( numFaces == size( topologyId, dim, 1 ) );

      delete[] origins;

      return numFaces;
    }

  } // namespace GenericGeometry

} // namespace Dune

#endif // DUNE_GEOMETRY_GENERICGEOMETRY_REFERENCEDOMAIN_HH
