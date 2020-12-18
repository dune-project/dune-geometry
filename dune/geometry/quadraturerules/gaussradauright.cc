// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include "config.h"
#include "../quadraturerules.hh"

namespace Dune {

  template GaussRadauRightQuadratureRule1D<float>::GaussRadauRightQuadratureRule1D(int);
  template GaussRadauRightQuadratureRule1D<double>::GaussRadauRightQuadratureRule1D(int);

} // namespace Dune
