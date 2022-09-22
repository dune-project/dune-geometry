// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#include <config.h>

#include <iostream>

#include <dune/geometry/type.hh>

int main ( int /* argc */, char ** /* argv */ )
{
  int fail = 0;
  for( int dim = 0; dim < 10; ++dim )
  {
    Dune::GeometryType gt = Dune::GeometryTypes::none( dim );

    if ( ! gt.isNone() )
    {
      fail = 1;
      std::cerr << "Geometry types 'none' for dim " << dim << " fails using makeNone " << std::endl;
    }
  }
  return fail;
}
