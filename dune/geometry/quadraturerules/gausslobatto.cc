// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include "config.h"
#include "../quadraturerules.hh"

namespace Dune {

  template GaussLobattoQuadratureRule1D<float>::GaussLobattoQuadratureRule1D(int);
  template GaussLobattoQuadratureRule1D<double>::GaussLobattoQuadratureRule1D(int);

} // namespace Dune
