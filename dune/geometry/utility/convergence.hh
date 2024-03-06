// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GEOMETRY_UTILITY_CONVERGENCE_HH
#define DUNE_GEOMETRY_UTILITY_CONVERGENCE_HH

#include <cmath>
#include <limits>

namespace Dune::Impl {

template <class R = double>
struct ConvergenceOptions
{
  //! Maximal number of iterations
  int maxIt = 100;

  //! Absolute tolerance for convergence
  R absTol = []{ using std::sqrt; return sqrt(std::numeric_limits<R>::epsilon()); }();
};

} // end namespace Dune::Impl

#endif // DUNE_GEOMETRY_UTILITY_CONVERGENCE_HH
