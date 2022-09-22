// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GEOMETRY_QUADRATURE_POINT_HH
#define DUNE_GEOMETRY_QUADRATURE_POINT_HH

#ifndef DUNE_INCLUDING_IMPLEMENTATION
#error This is a private header that should not be included directly.
#error Use #include <dune/geometry/quadraturerules.hh> instead.
#endif

namespace Dune {

  /** \brief Quadrature for a point (0D) */
  template<typename ct>
  class PointQuadratureRule :
    public QuadratureRule<ct,0>
  {
    // compile time parameters
    constexpr static int dim = 0;

    friend class QuadratureRuleFactory<ct,dim>;

    PointQuadratureRule () : QuadratureRule<ct,0>(GeometryTypes::vertex)
    {
      FieldVector<ct, dim> point(0.0);

      // Any function is integrated exactly, hence the order is infinite
      this->delivered_order = std::numeric_limits<int>::max();
      this->push_back(QuadraturePoint<ct,dim>(point, 1.0));
    }

    ~PointQuadratureRule(){}

  };

} // end namespace Dune

#endif // DUNE_GEOMETRY_QUADRATURE_POINT_HH
