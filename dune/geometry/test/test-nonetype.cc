#include <config.h>

#include <iostream>

#include <dune/geometry/type.hh>

int main ( int argc, char **argv )
{
  int fail = 0;
  for( int dim = 0; dim < 10; ++dim )
  {
    Dune::GeometryType gt1( Dune::GeometryType::none, dim );

    Dune::GeometryType gt2;
    gt2.makeNone( dim );

    if( gt1 != gt2 )
    {
      fail = 1;
      std::cerr << "Geometry types 'none' for dim " << dim << " do not coincide" << std::endl;
    }
  }
  return fail;
};
