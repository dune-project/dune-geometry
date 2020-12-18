// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include "config.h"
#include "../quadraturerules.hh"

namespace Dune {

  template GaussRadauLeftQuadratureRule1D<float>::GaussRadauQuadratureRule1D(int);
  template GaussRadauLeftQuadratureRule1D<double>::GaussRadauQuadratureRule1D(int);

} // namespace Dune
