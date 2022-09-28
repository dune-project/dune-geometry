// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
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
  template class GaussLobattoQuadratureRule<double, 1>;
  template class GaussQuadratureRule<double, 1>;
  template class GaussRadauLeftQuadratureRule<double, 1>;
  template class GaussRadauRightQuadratureRule<double, 1>;
  template class Jacobi1QuadratureRule<double, 1>;
  template class Jacobi2QuadratureRule<double, 1>;
  template class JacobiNQuadratureRule<double, 1>;
  template class PrismQuadratureRule<double, 3>;
  template class SimplexQuadratureRule<double, 2>;
  template class SimplexQuadratureRule<double, 3>;

} // namespace
