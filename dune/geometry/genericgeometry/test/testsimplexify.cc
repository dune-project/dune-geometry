// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>
#include <vector>

#include <dune/geometry/genericgeometry/simplexify.hh>

const unsigned int tetrahedron[ 1 ][ 4 ] = { { 0, 1, 2, 3 } };
const unsigned int pyramid[ 2 ][ 4 ] = { { 0, 1, 3, 4 }, { 0, 2, 3, 4 } };
const unsigned int prism[ 3 ][ 4 ] = { { 0, 1, 2, 5 }, { 0, 1, 4, 5 }, { 0, 3, 4, 5 } };
const unsigned int hexahedron[ 6 ][ 4 ] = { { 0, 1, 3, 7 }, { 0, 2, 3, 7 }, { 0, 1, 5, 7 },
                                            { 0, 2, 6, 7 }, { 0, 4, 5, 7 }, { 0, 4, 6, 7 } };

bool test ( unsigned int topologyId, unsigned int numSimplices, const unsigned int (*correct)[ 4 ] )
{
  std::vector< std::vector< unsigned int > > simplices;
  Dune::GenericGeometry::simplexify( topologyId, 3, simplices );

  if( simplices.size() != numSimplices )
  {
    std::cerr << "Error for topologyId " << topologyId << ": Wrong number of simplices ";
    std::cerr << "(" << simplices.size() << ", should be " << numSimplices << ")." << std::endl;
    return false;
  }

  for( unsigned int i = 0; i < numSimplices; ++i )
  {
    if( simplices[ i ].size() != 4 )
    {
      std::cerr << "Error for topologyId " << topologyId << ": Wrong number of vertices for simplex ";
      std::cerr << "(" << simplices[ i ].size() << ")." << std::endl;
      return false;
    }

    for( int j = 0; j <= 3; ++j )
    {
      if( simplices[ i ][ j ] != correct[ i ][ j ] )
      {
        std::cerr << "Error for topologyId " << topologyId << ": Wrong simplex returned ";
        std::cerr << "( [" << simplices[ i ][ 0 ];
        for( int k = 1; k <= 3; ++k )
          std::cerr << ", " << simplices[ i ][ k ];
        std::cerr << "], should be [" << correct[ i ][ 0 ];
        for( int k = 1; k <= 3; ++k )
          std::cerr << ", " << correct[ i ][ k ];
        std::cerr <<  "])." << std::endl;
        return false;
      }
    }

  }
  return true;
}

int main ( int argc, char **argv )
{
  bool success = true;
  success &= test( 0, 1, tetrahedron );
  success &= test( 1, 1, tetrahedron );
  success &= test( 2, 2, pyramid );
  success &= test( 3, 2, pyramid );
  success &= test( 4, 3, prism );
  success &= test( 5, 3, prism );
  success &= test( 6, 6, hexahedron );
  success &= test( 7, 6, hexahedron );

  return (success ? 0 : 1);
}
