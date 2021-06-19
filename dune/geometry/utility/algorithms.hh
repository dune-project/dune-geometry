// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GEOMETRY_UTILITY_ALGORITHMS_HH
#define DUNE_GEOMETRY_UTILITY_ALGORITHMS_HH

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <type_traits>

#include <dune/common/debugstream.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/ftraits.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/affinegeometry.hh> // for FieldMatrixHelper

namespace Dune {
namespace Impl {

template <class R = double>
struct GaussNewtonOptions
{
  //! Maximal number of Newton iterations
  int maxIt = 100;

  //! Break tolerance for the residual |f(x) - y| in the 2-norm.
  R absTol = []{ using std::sqrt; return sqrt(std::numeric_limits<R>::epsilon()); }();

  //! Maximal number of line-search iterations
  int maxInnerIt = 10;

  //! Reduction factor for the line-search parameter
  R theta = 0.5;
};


/// \brief A return-code indicating the outcome of the `gaussNewton` algorithm.
enum class GaussNewtonErrorCode
{
  OK = 0,                   //< A solution is found
  JACOBIAN_NOT_INVERTIBLE,  //< The Jacobian is not invertible at the current point
  STAGNATION,               //< No reduction of the residul norm possible
  TOLERANCE_NOT_REACHED     //< The break tolerance for the resodual norm is not reached
};


/**
 * \brief Gauß-Newton method to solve the least-squares problem `|f(x) - y|^2 -> min`
 *
 * \param[in] f     Objective function
 * \param[in] df    Gradient of the objective function f
 * \param[in] y     Target value
 * \param[inout] x0 Initial guess for the solution and the result
 * \param[in] opts  Additional parameter to control the method, see \ref `GaussNewtonOptions`
 * \result Argument `x0` such that `f(x0) = y` in the sense of least squares error.
 *         If a solution could be found the return code `GaussNewtonErrorCode::OK`
 *         is returned, otherwise an indicator for the error.
 */
template <class F, class DF, class Domain,
          class Range = std::invoke_result_t<F, Domain>,
          class R = typename Dune::FieldTraits<Domain>::real_type>
GaussNewtonErrorCode gaussNewton (const F& f, const DF& df, Range y, Domain& x0,
                                  GaussNewtonOptions<R> opts = {})
{
  Domain x = x0;
  Domain dx{};
  Range dy = f(x0) - y;
  R resNorm0 = dy.two_norm();
  R resNorm = 0;

  for (int i = 0; i < opts.maxIt; ++i)
  {
    // Get descent direction dx: (J^T*J)dx = J^T*dy
    const bool invertible = FieldMatrixHelper<R>::xTRightInvA(df(x), dy, dx);

    // break if jacobian is not invertible
    if (!invertible)
      return GaussNewtonErrorCode::JACOBIAN_NOT_INVERTIBLE;

    // line-search procedure to update x with correction dx
    R alpha = 1;
    for (int j = 0; j < opts.maxInnerIt; ++j) {
      x = x0 - alpha * dx;
      dy = f(x) - y;
      resNorm = dy.two_norm();

      if (resNorm < resNorm0)
        break;

      alpha *= opts.theta;
    }

    // cannot reduce the residual
    if (!(resNorm < resNorm0))
      return GaussNewtonErrorCode::STAGNATION;

    x0 = x;
    resNorm0 = resNorm;

    // break if tolerance is reached.
    if (resNorm < opts.absTol)
      return GaussNewtonErrorCode::OK;
  }

  // tolerance could not be reached
  if (!(resNorm < opts.absTol))
    return GaussNewtonErrorCode::TOLERANCE_NOT_REACHED;

  return GaussNewtonErrorCode::OK;
}

} // end namespace Impl
} // end namespace Dune

#endif // DUNE_GEOMETRY_UTILITY_ALGORITHMS_HH
