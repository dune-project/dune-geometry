// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_GENERICGEOMETRY_MATRIXHELPER_HH
#define DUNE_GEOMETRY_GENERICGEOMETRY_MATRIXHELPER_HH

#warning This header is deprecated. Please #include <dune/common/matrixhelper.hh> instead

#include <dune/common/matrixhelper.hh>

namespace Dune
{

  namespace GenericGeometry
  {
    template< class Field >
    using FieldHelper = Dune::FieldHelper<Field>;

    template< class Traits >
    using MatrixHelper = Dune::MatrixHelper<Traits>;
  } // namespace GenericGeometry

} // namespace Dune

#endif // #ifndef DUNE_GEOMETRY_GENERICGEOMETRY_MATRIXHELPER_HH
