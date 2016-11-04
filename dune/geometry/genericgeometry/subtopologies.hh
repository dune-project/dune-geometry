// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_GENERICGEOMETRY_SUBTOPOLOGIES_HH
#define DUNE_GEOMETRY_GENERICGEOMETRY_SUBTOPOLOGIES_HH

#include <cassert>
#include <vector>

#include <dune/common/forloop.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/visibility.hh>
#include <dune/common/unused.hh>

#include <dune/geometry/genericgeometry/topologytypes.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    /** \brief Compute the number of subentities of a given codimension */
    unsigned int size ( unsigned int topologyId, int dim, int codim );



    /** \brief Compute the topology id of a given subentity
     *
     * \param topologyId Topology id of entity
     * \param dim Dimension of entity
     * \param codim Codimension of the subentity that we are interested in
     * \param i Number of the subentity that we are interested in
     */
    unsigned int subTopologyId ( unsigned int topologyId, int dim, int codim, unsigned int i );



    // subTopologyNumbering
    // --------------------

    void subTopologyNumbering ( unsigned int topologyId, int dim, int codim, unsigned int i, int subcodim,
                                unsigned int *beginOut, unsigned int *endOut );

  } // namespace GenericGeometry

} // namespace Dune

#endif // #ifndef DUNE_GEOMETRY_GENERICGEOMETRY_SUBTOPOLOGIES_HH
