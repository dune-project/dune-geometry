// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#include "config.h"
#include "../quadraturerules.hh"

namespace Dune {

  template Jacobi2QuadratureRule1D<float>::Jacobi2QuadratureRule1D(int);
  template Jacobi2QuadratureRule1D<double>::Jacobi2QuadratureRule1D(int);

} // namespace
