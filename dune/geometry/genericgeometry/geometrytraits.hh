// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_GENERICGEOMETRY_GEOMETRYTRAITS_HH
#define DUNE_GEOMETRY_GENERICGEOMETRY_GEOMETRYTRAITS_HH

#include "../type.hh"
#include "matrixhelper.hh"

namespace Dune
{
  namespace GenericGeometry
  {

    // DuneCoordTraits
    // ---------------

    template< class ct >
    struct DuneCoordTraits
    {
      typedef ct ctype;

      template< int dim >
      struct Vector
      {
        typedef FieldVector< ctype, dim > type;
      };

      template< int rows, int cols >
      struct Matrix
      {
        typedef FieldMatrix< ctype, rows, cols > type;
      };

      // This limit is, e.g., used in the termination criterion of the Newton
      // scheme within the generic implementation of the method local
      static const ctype epsilon ()
      {
        return 1e-6;
      }
    };

  } // namespace GenericGeometry

} // namespace Dune

#endif // #ifndef DUNE_GEOMETRY_GENERICGEOMETRY_GEOMETRYTRAITS_HH
