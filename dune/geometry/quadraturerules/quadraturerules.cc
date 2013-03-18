// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"
#include "../quadraturerules.hh"

namespace Dune {

  /** Singleton holding the Gauss points on the interval */
  SimplexQuadraturePoints<2> SimplexQuadraturePointsSingleton<2>::sqp;

  /** Singleton holding the SimplexQuadrature points dim==3 */
  SimplexQuadraturePoints<3> SimplexQuadraturePointsSingleton<3>::sqp;

  /** Singleton holding the Prism Quadrature points  */
  PrismQuadraturePoints<3> PrismQuadraturePointsSingleton<3>::prqp;

  template SimplexQuadratureRule<float, 2>::SimplexQuadratureRule(int);
  template SimplexQuadratureRule<double, 2>::SimplexQuadratureRule(int);
  template SimplexQuadratureRule<float, 3>::SimplexQuadratureRule(int);
  template SimplexQuadratureRule<double, 3>::SimplexQuadratureRule(int);

} // namespace
