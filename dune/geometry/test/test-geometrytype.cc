// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <string>
#include <iostream>
#include <vector>

#include <dune/geometry/type.hh>
#include <dune/geometry/genericgeometry/topologytypes.hh>
#include <dune/geometry/genericgeometry/subtopologies.hh>

std::string convBase(unsigned long v, long base)
{
  const char* digits = "0123456789abcdef";
  std::string result;
  if((base < 2) || (base > 16)) {
    result = "Error: base out of range.";
  }
  else {
    do {
      result = digits[v % base] + result;
      v /= base;
    }
    while(v);
  }
  return result;
}

unsigned int guessTopologyId(unsigned int dim, unsigned int vertices)
{
  if( (vertices < dim+1) || (vertices > (1u << dim)) )
    DUNE_THROW( Dune::Exception, "IMPOSSIBLE" );
  std::vector< unsigned int > ids;
  for( unsigned int id = 0; id < Dune::GenericGeometry::numTopologies( dim ); id += 2u )
  {
    if( Dune::GenericGeometry::size( id, dim, dim ) == vertices )
      ids.push_back( id );
  }
  if( ids.empty() )
    DUNE_THROW( Dune::Exception, "Impossible setting" );
  for( size_t i = 0; i < ids.size(); ++i )
    std::cout << "possibility " << i << "\t" << convBase( ids[ i ], 2 ) << std::endl;
  if( ids.size() > 1 )
    DUNE_THROW( Dune::Exception, "Too many options" );
  return ids[ 0 ];
}

void testGuess(unsigned int dim, unsigned int vertices)
{
  std::cout << "dim: " << dim << std::endl;
  std::cout << "vertices: " << vertices << std::endl;
  unsigned int id = guessTopologyId(dim, vertices);
  std::cout << "guess: " << convBase(id, 2) << std::endl;
  // Dune::GeometryType gt;
  // gt.makeFromVertices(dim, vertices);
  // std::cout << "real:  " << convBase(gt.id(), 2) << std::endl;
}

int main()
{
  for (int d=0; d<=8; d++)
    for (int v=d+1; v<=(1<<d); v++)
    {
      try {
        testGuess(d,v);
      }
      catch (Dune::Exception & e)
      {
        std::cout << "Error: " << e.what() << std::endl;
      }
    }
}
