// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
//

#ifndef DUNE_INCLUDING_IMPLEMENTATION
#error This is a private header that should not be included directly.
#error Use #include <dune/geometry/quadraturerules.hh> instead.
#endif
#undef DUNE_INCLUDING_IMPLEMENTATION

namespace Dune {

  // for fundamental types
  template<typename ct>
  void GaussRadauLeftQuadratureInitHelper<ct,true>::init(int p,
         std::vector< FieldVector<ct, 1> > & _points,
         std::vector< ct > & _weight,
         int & delivered_order)
  {
    switch(p)
    {
    // order 0,1
    case 0 :
        delivered_order = 0.;
        _points.resize(1);
        _weight.resize(1);
        _points[0] = 0. ;
        _weight[0] = 1. ;
        break;
    case 1 :
      delivered_order = 1;
      _points.resize(2);
      _weight.resize(2);
      _points[0] = 0.0000000000000000e+00;
      _weight[0] = 2.5000000000000000e-01;
      _points[1] = 6.6666666666666663e-01;
      _weight[1] = 7.5000000000000000e-01;
      break;

    // order 2,3
    case 2 :
    case 3 :
      delivered_order = 3;
      _points.resize(3);
      _weight.resize(3);
      _points[0] = 0.0;
      _weight[0] = 1.1111111111111110e-01;
      _points[1] = 3.5505102572168223e-01;
      _weight[1] = 5.1248582618842164e-01;
      _points[2] = 8.4494897427831783e-01;
      _weight[2] = 3.7640306270046725e-01;
      break;
    default :
      DUNE_THROW(QuadratureOrderOutOfRange, "Quadrature rule " << p << " not supported!");
    }
  }

  // for non-fundamental types: assign numbers as strings
  template<typename ct>
  void GaussRadauLeftQuadratureInitHelper<ct,false>::init(int p,
         std::vector< FieldVector<ct, 1> > & _points,
         std::vector< ct > & _weight,
         int & delivered_order)
  {
    switch(p)
    {
    default :
      DUNE_THROW(QuadratureOrderOutOfRange, "Quadrature rule " << p << " not supported!");
    }
  }

} // namespace
