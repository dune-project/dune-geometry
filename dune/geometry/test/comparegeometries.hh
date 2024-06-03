// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GEOMETRY_TEST_COMPAREGEOMETRIES_HH
#define DUNE_GEOMETRY_TEST_COMPAREGEOMETRIES_HH

#include <cmath>
#include <iostream>
#include <limits>
#include <type_traits>

namespace Dune {

template <class R>
R defaultTolerance ()
{
  using std::sqrt;
  return sqrt(std::numeric_limits<R>::epsilon());
}

template <class Geometry1, class Geometry2,
  class R = std::common_type_t<typename Geometry1::ctype, typename Geometry2::ctype>>
bool compareGeometries (const Geometry1& geo1, const Geometry2& geo2,
                        R tolerance = defaultTolerance<R>())
{
  if constexpr(Geometry1::mydimension != Geometry2::mydimension) {
    std::cerr << "Error: Dimensions do not match." << std::endl;
    return false;
  }
  else if constexpr(Geometry1::coorddimension != Geometry2::coorddimension) {
    std::cerr << "Error: Coord-dimensions do not match." << std::endl;
    return false;
  }
  else {
    using std::abs;

    if (geo1.affine() != geo2.affine()) {
      std::cerr << "Error: Affine-property does not match." << std::endl;
      return false;
    }

    if (geo1.type() != geo2.type()) {
      std::cerr << "Error: GeometryType does not match." << std::endl;
      return false;
    }

    if (geo1.corners() != geo2.corners()) {
      std::cerr << "Error: Number of corners does not match." << std::endl;
      return false;
    }

    for (int i = 0; i < geo1.corners(); ++i) {
      if ((geo1.corner(i) - geo2.corner(i)).two_norm() > tolerance) {
        std::cerr << "Error: Corner " << i << " does not match." << std::endl;
        return false;
      }
    }

    if ((geo1.center() - geo2.center()).two_norm() > tolerance) {
      std::cerr << "Error: Center does not match." << std::endl;
      return false;
    }

    if (abs(geo1.volume() - geo2.volume()) > tolerance) {
      std::cerr << "Error: Volume does not match." << std::endl;
      return false;
    }

    const auto& quadrature = Dune::QuadratureRules<R, Geometry1::mydimension>::rule(geo1.type(), 4);
    for (auto&& [pos,weight] : quadrature)
    {
      if ((geo1.global(pos) - geo2.global(pos)).two_norm() > tolerance) {
        std::cerr << "Error: global(" << pos << ") does not match." << std::endl;
        return false;
      }

      if ((geo1.jacobian(pos) - geo2.jacobian(pos)).frobenius_norm() > tolerance) {
        std::cerr << "Error: jacobian(" << pos << ") does not match." << std::endl;
        return false;
      }

      if ((geo1.jacobianTransposed(pos) - geo2.jacobianTransposed(pos)).frobenius_norm() > tolerance) {
        std::cerr << "Error: jacobianTransposed(" << pos << ") does not match." << std::endl;
        return false;
      }

      if ((geo1.jacobianInverse(pos) - geo2.jacobianInverse(pos)).frobenius_norm() > tolerance) {
        std::cerr << "Error: jacobianInverse(" << pos << ") does not match." << std::endl;
        return false;
      }

      if ((geo1.jacobianInverseTransposed(pos) - geo2.jacobianInverseTransposed(pos)).frobenius_norm() > tolerance) {
        std::cerr << "Error: jacobianInverse(" << pos << ") does not match." << std::endl;
        return false;
      }

      if (abs(geo1.integrationElement(pos) - geo2.integrationElement(pos)) > tolerance) {
        std::cerr << "Error: integrationElement(" << pos << ") does not match." << std::endl;
        return false;
      }
    }

    return true;
  }
}

} // end namespace Dune

#endif // DUNE_GEOMETRY_TEST_COMPAREGEOMETRIES_HH
