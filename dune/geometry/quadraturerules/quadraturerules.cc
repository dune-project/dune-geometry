// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"

#define DUNE_NO_EXTERN_QUADRATURERULES
#include "../quadraturerules.hh"
#undef DUNE_NO_EXTERN_QUADRATURERULES

namespace Dune {

  /** Singleton holding the Gauss points on the interval */
  SimplexQuadraturePoints<2> SimplexQuadraturePointsSingleton<2>::sqp;

  /** Singleton holding the SimplexQuadrature points dim==3 */
  SimplexQuadraturePoints<3> SimplexQuadraturePointsSingleton<3>::sqp;

  /** Singleton holding the Prism Quadrature points  */
  PrismQuadraturePoints<3> PrismQuadraturePointsSingleton<3>::prqp;

  // explicit template instatiation
  template class Jacobi1QuadratureRule<float, 1>;
  template class Jacobi1QuadratureRule<double, 1>;
  template class Jacobi2QuadratureRule<float, 1>;
  template class Jacobi2QuadratureRule<double, 1>;
  template class GaussQuadratureRule<float, 1>;
  template class GaussQuadratureRule<double, 1>;
  template class GaussLobattoQuadratureRule<float, 1>;
  template class GaussLobattoQuadratureRule<double, 1>;

  template class SimplexQuadratureRule<float, 2>;
  template class SimplexQuadratureRule<double, 2>;
  template class SimplexQuadratureRule<float, 3>;
  template class SimplexQuadratureRule<double, 3>;

  template class PrismQuadratureRule<float, 3>;
  template class PrismQuadratureRule<double, 3>;

} // namespace
