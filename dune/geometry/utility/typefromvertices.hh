// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_TYPE_FROM_VERTICES_HH
#define DUNE_GEOMETRY_TYPE_FROM_VERTICES_HH

#include <dune/geometry/type.hh>

namespace Dune {

  /** \brief Utitlity function to construct the correct geometry type given the dimension and the number of vertices
   *  \note This code only works up to dimension 3.
   *        In higher dimensions the number of vertices does not uniquely identify the type of polyhedron.
   */
  GeometryType geometryTypeFromVertices(unsigned int dim, unsigned int vertices)
  {
    GeometryType gt;
    switch (dim)
    {
      case 0 :
        gt.makeVertex();
        return gt;
      case 1 :
        gt.makeLine();
        return gt;
      case 2 :
        switch (vertices) {
          case 3 :
            gt.makeSimplex(2);
            return gt;
          case 4 :
            gt.makeCube(2);
            return gt;
          default :
            DUNE_THROW(NotImplemented, "2d elements with " << vertices << " corners are not supported!");
        }
      case 3 :
        switch (vertices) {
          case 4 :
            gt.makeSimplex(3);
            return gt;
          case 5 :
            gt.makePyramid();
            return gt;
          case 6 :
            gt.makePrism();
            return gt;
          case 8 :
            gt.makeCube(3);
            return gt;
          default :
            DUNE_THROW(NotImplemented, "3d elements with " << vertices << " corners are not supported!");
        }
      default :
        DUNE_THROW(NotImplemented, "geometryTypeFromVertices works only up to dim=3");
    }
  }

}

#endif // DUNE_GEOMETRY_TYPE_FROM_VERTICES_HH
