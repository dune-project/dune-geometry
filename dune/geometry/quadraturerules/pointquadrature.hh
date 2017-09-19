// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
namespace Dune {

  /** \brief Quadrature for a point (0D) */
  template<typename ct>
  class PointQuadratureRule :
    public QuadratureRule<ct,0>
  {
    // compile time parameters
    enum { dim=0 };

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
