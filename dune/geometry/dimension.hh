// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_DIMENSION_HH
#define DUNE_GEOMETRY_DIMENSION_HH

#include <type_traits>

namespace Dune {

  //! Static tag representing a dimension.
  template<int dim>
  struct Dim
    : public std::integral_constant<int,dim>
  {
    typedef Dim type;
  };

  //! Static tag representing a codimension.
  template<int codim>
  struct Codim
    : public std::integral_constant<int,codim>
  {
    typedef Codim type;
  };

}

#endif // DUNE_GEOMETRY_DIMENSION_HH
