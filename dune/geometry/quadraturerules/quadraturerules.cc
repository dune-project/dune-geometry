// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"

#define DUNE_NO_EXTERN_QUADRATURERULES
#include "../quadraturerules.hh"
#undef DUNE_NO_EXTERN_QUADRATURERULES

namespace Dune {

  // explicit template instatiation
  template class Jacobi1QuadratureRule<float, 1>;
  template class Jacobi1QuadratureRule<double, 1>;
  template class Jacobi2QuadratureRule<float, 1>;
  template class Jacobi2QuadratureRule<double, 1>;
  template class GaussQuadratureRule<float, 1>;
  template class GaussQuadratureRule<double, 1>;
  template class GausslobattoQuadratureRule<float, 1>;
  template class GausslobattoQuadratureRule<double, 1>;

  template class Cube2dQuadratureRule<float, 2>;
  template class Cube2dQuadratureRule<double, 2>;
  template class Simplex2dQuadratureRule<float, 2>;
  template class Simplex2dQuadratureRule<double, 2>;

  template class Cube3dQuadratureRule<float, 3>;
  template class Cube3dQuadratureRule<double, 3>;
  template class Simplex3dQuadratureRule<float, 3>;
  template class Simplex3dQuadratureRule<double, 3>;
  template class PrismQuadratureRule<float, 3>;
  template class PrismQuadratureRule<double, 3>;
  template class PyramidQuadratureRule<float, 3>;
  template class PyramidQuadratureRule<double, 3>;

} // namespace
