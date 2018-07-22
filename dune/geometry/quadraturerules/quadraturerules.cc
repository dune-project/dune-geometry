// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"

#define DUNE_NO_EXTERN_QUADRATURERULES
#include "../quadraturerules.hh"
#undef DUNE_NO_EXTERN_QUADRATURERULES

namespace Dune {

  // explicit template instatiation
  template class Jacobi1QuadratureRule<double>;
  template class Jacobi2QuadratureRule<double>;
  template class GaussQuadratureRule<double>;
  template class GaussLobattoQuadratureRule<double>;
  template class PrismQuadratureRule<double>;
  template class SimplexQuadratureRule<double, 2>;
  template class SimplexQuadratureRule<double, 3>;

} // namespace
