// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include "config.h"
#include "../quadraturerules.hh"

namespace Dune {

  template Jacobi1QuadratureRule1D<float>::Jacobi1QuadratureRule1D(int);
  template Jacobi1QuadratureRule1D<double>::Jacobi1QuadratureRule1D(int);

} // namespace
