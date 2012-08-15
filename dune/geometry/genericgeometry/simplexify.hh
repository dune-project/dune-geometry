// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_GENERICGEOMETRY_SIMPLEXIFY_HH
#define DUNE_GEOMETRY_GENERICGEOMETRY_SIMPLEXIFY_HH

#include <vector>

namespace Dune
{

  namespace GenericGeometry
  {

    void simplexify ( unsigned int topologyId, int dim, std::vector< std::vector< unsigned int > > &simplices );

  } // namespace GenericGeometry

} // namespace Dune

#endif // #ifndef DUNE_GEOMETRY_GENERICGEOMETRY_SIMPLEXIFY_HH
