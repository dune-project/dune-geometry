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

    if ( ! gt1.isNone() )
    {
      fail = 1;
      std::cerr << "Geometry types 'none' for dim " << dim << " has wrong constructor " << std::endl;
    }

    if ( ! gt2.isNone() )
    {
      fail = 1;
      std::cerr << "Geometry types 'none' for dim " << dim << " fails using makeNone " << std::endl;
    }

    if( gt1 != gt2 )
    {
      fail = 1;
      std::cerr << "Geometry types 'none' for dim " << dim << " do not coincide" << std::endl;
    }
  }
  return fail;
}
