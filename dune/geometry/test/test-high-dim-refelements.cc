// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/geometry/referenceelements.hh>

typedef void ctype;
static const int dim = 5;

int main ( int argc, char **argv )
{
  const Dune::GeometryType gt( Dune::GeometryType::cube, dim );
  const Dune::ReferenceElement< ctype, dim > &refElement
    = Dune::ReferenceElements< ctype, dim >::general( gt );

  std::cout << "number of vertices: " << refElement.size( dim ) << std::endl;

  return 0;
}
