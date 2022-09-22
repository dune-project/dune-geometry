// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GEOMETRY_TYPE_FROM_VERTEX_COUNT_HH
#define DUNE_GEOMETRY_TYPE_FROM_VERTEX_COUNT_HH

#include <dune/geometry/type.hh>

namespace Dune {

  /** \brief Utitlity function to construct the correct geometry type given the dimension and the number of vertices
   *  \note This code only works up to dimension 3.
   *        In higher dimensions the number of vertices does not uniquely identify the type of polyhedron.
   */
  inline
  GeometryType geometryTypeFromVertexCount(unsigned int dim, unsigned int vertices)
  {
    switch (dim)
    {
      case 0 :
        return GeometryTypes::vertex;
      case 1 :
        return GeometryTypes::line;
      case 2 :
        switch (vertices) {
          case 3 :
            return GeometryTypes::triangle;
          case 4 :
            return GeometryTypes::quadrilateral;
          default :
            DUNE_THROW(NotImplemented, "2d elements with " << vertices << " corners are not supported!");
        }
      case 3 :
        switch (vertices) {
          case 4 :
            return GeometryTypes::tetrahedron;
          case 5 :
            return GeometryTypes::pyramid;
          case 6 :
            return GeometryTypes::prism;
          case 8 :
            return GeometryTypes::hexahedron;
          default :
            DUNE_THROW(NotImplemented, "3d elements with " << vertices << " corners are not supported!");
        }
      default :
        DUNE_THROW(NotImplemented, "geometryTypeFromVertexCount works only up to dim=3");
    }
  }

}

#endif // DUNE_GEOMETRY_TYPE_FROM_VERTEX_COUNT_HH
