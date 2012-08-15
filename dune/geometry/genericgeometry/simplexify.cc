// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/geometry/genericgeometry/simplexify.hh>
#include <dune/geometry/genericgeometry/topologytypes.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // numSimplices
    // ------------

    static std::size_t numSimplices ( unsigned int topologyId, int dim )
    {
      assert( (dim >= 0) && (topologyId < numTopologies( dim )) );

      if( dim > 0 )
      {
        if( isPrism( topologyId, dim ) )
          return dim * numSimplices( baseTopologyId( topologyId, dim ), dim-1 );
        else
          return numSimplices( baseTopologyId( topologyId, dim ), dim-1 );
      }
      else
        return 1;
    }



    // doSimplexify
    // ------------

    static std::size_t doSimplexify ( unsigned int topologyId, int dim, std::vector< std::vector< unsigned int > > &simplices )
    {
      assert( (dim >= 0) && (topologyId < numTopologies( dim )) );

      if( dim > 0 )
      {
        std::size_t baseSize = doSimplexify( baseTopologyId( topologyId, dim ), dim-1, simplices );
        // note: the last vertex of all simplices will have the "largest"
        const unsigned int numBaseVertices = simplices[ 0 ][ dim-1 ] + 1;

        if( isPrism( topologyId, dim ) )
        {
          assert( simplices.size() >= dim*baseSize );
          for( int k = 0; k < dim; ++k )
          {
            for( std::size_t i = 0; i < baseSize; ++i )
            {
              const std::size_t idx = k*baseSize + i;
              assert( idx < simplices.size() );
              for( int j = 0; j <= dim; ++j )
              {
                const int s = (j < dim-k) ? 0 : 1;
                simplices[ idx ][ j ] = s*numBaseVertices + simplices[ i ][ j-s ];
              }
            }
          }
          return dim * baseSize;
        }
        else
        {
          assert( isPyramid( topologyId, dim ) );
          for( std::size_t i = 0; i < baseSize; ++i )
          {
            assert( simplices[ i ].size() > std::size_t( dim ) );
            simplices[ i ][ dim ] = numBaseVertices;
          }
          return baseSize;
        }
      }
      else
      {
        assert( simplices[ 0 ].size() > std::size_t( dim ) );
        simplices[ 0 ][ 0 ] = 0;
        return 1;
      }
    }



    // simplexify
    // ----------

    void simplexify ( unsigned int topologyId, int dim, std::vector< std::vector< unsigned int > > &simplices )
    {
      simplices.resize( numSimplices( topologyId, dim ), std::vector< unsigned int >( dim+1 ) );
      doSimplexify( topologyId, dim, simplices );
    }

  } // namespace GenericGeometry

} // namespace Dune
