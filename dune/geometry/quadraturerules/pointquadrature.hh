// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_QUADRATURE_POINT_HH
#define DUNE_GEOMETRY_QUADRATURE_POINT_HH

#ifndef DUNE_INCLUDING_IMPLEMENTATION
#error This is a private header that should not be included directly.
#error Use #include <dune/geometry/quadraturerules.hh> instead.
#endif

namespace Dune {

  /** \brief Quadrature for a point (0D) */
  template<typename ct>
  class PointQuadratureRule : public QuadratureRule<ct,0>
  {
    friend class QuadratureRuleFactory<ct,0>;

    PointQuadratureRule ()
      : QuadratureRule<ct,0>(GeometryTypes::vertex, std::numeric_limits<int>::max())
    {
      this->emplace_back(FieldVector<ct,0>{}, 1.0);
    }
  };

} // end namespace Dune

#endif // DUNE_GEOMETRY_QUADRATURE_POINT_HH
