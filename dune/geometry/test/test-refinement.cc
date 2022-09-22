// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

/*!
 * \file
 * \brief Unit tests for the virtual refinement code
 */

#include "config.h"

#include <iostream>
#include <ostream>
#include <typeinfo>

#include <dune/geometry/test/checkgeometry.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/virtualrefinement.hh>

using namespace Dune;

void pass(int &result)
{
  // 77 means 'SKIP'
  if (result == 77)
  {
    result = 0;
  }
}

void fail(int &result)
{
  result = 1;
}

void collect(int &result, bool passed)
{
  if (passed)
  {
    pass(result);
  }
  else
  {
    fail(result);
  }
}

/*!
 * \brief Test virtual refinement for an element with a run-time type
 */
template <class ct, int dim>
void testVirtualRefinement(int &result, const Dune::GeometryType& elementType,
                           const Dune::GeometryType& coerceTo, Dune::RefinementIntervals tag, std::string refType)
{
  std::cout << "Checking virtual refinement " << elementType << " -> "
            << coerceTo << " intervals " << tag.intervals()
            << " " << refType << " tag" << std::endl;

  auto refElem = referenceElement<ct, dim>(elementType);

  typedef Dune::VirtualRefinement<dim, ct> Refinement;
  typedef typename Refinement::ElementIterator eIterator;
  typedef typename Refinement::VertexIterator vIterator;

  // Make a virtual refinement of the reference element
  Refinement & elementRefinement =
    Dune::buildRefinement<dim, ct>(elementType, coerceTo);

  eIterator eSubEnd = elementRefinement.eEnd(tag);
  eIterator eSubIt  = elementRefinement.eBegin(tag);

  for (; eSubIt != eSubEnd; ++eSubIt)
  {
    if (refElem.checkInside(eSubIt.coords()))
    {
      pass(result);
    }
    else
    {
      std::cerr << "Error: Sub-element position (" << eSubIt.coords()
                << ") is outside of the reference element" << std::endl;
      fail(result);
    }
  }

  vIterator vSubEnd = elementRefinement.vEnd(tag);
  vIterator vSubIt  = elementRefinement.vBegin(tag);

  for (; vSubIt != vSubEnd; ++vSubIt)
  {
    if (refElem.checkInside(vSubIt.coords()))
    {
      pass(result);
    }
    else
    {
      std::cerr << "Error: Sub-vertex position (" << vSubIt.coords()
                << ") is outside of the reference element" << std::endl;
      fail(result);
    }
  }
}

/*!
 * \brief Test virtual refinement for an element with a static type
 */
template <unsigned topologyId, class ct, unsigned coerceToId, int dim>
void testStaticRefinementGeometry(int &result, Dune::RefinementIntervals tag, std::string refType)
{
  std::cout << "Checking static refinement geometry "
            << GeometryType(topologyId, dim) << " -> "
            << GeometryType(coerceToId, dim) << " intervals " << tag.intervals()
            << " " << refType << " tag" << std::endl;

  typedef Dune::StaticRefinement<topologyId, ct, coerceToId, dim> Refinement;
  typedef typename Refinement::ElementIterator eIterator;
  typedef typename Refinement::VertexIterator vIterator;

  eIterator eSubEnd = Refinement::eEnd(tag);
  eIterator eSubIt  = Refinement::eBegin(tag);

  for (; eSubIt != eSubEnd; ++eSubIt)
  {
    // Call the standard test for geometries
    collect(result, checkGeometry(eSubIt.geometry()));
  }

  vIterator vSubEnd = Refinement::vEnd(tag);
  vIterator vSubIt  = Refinement::vBegin(tag);

  for (; vSubIt != vSubEnd; ++vSubIt)
  {
    // Call the standard test for geometries
    collect(result, checkGeometry(vSubIt.geometry()));
  }
}


int main(int /* argc */, char** /* argv */) try
{
  // 77 means 'SKIP'
  int result = 77;

  GeometryType gt1, gt2;

  // test segment
  gt1 = gt2 = GeometryTypes::line;
  for (unsigned int refinement = 0; refinement < 3; refinement++)
  {
    testVirtualRefinement<double,1>(result, gt1, gt2, refinementLevels(refinement), "levels");
    testVirtualRefinement<double,1>(result, gt1, gt2, refinementIntervals(1<<refinement), "intervals");
    testStaticRefinementGeometry<GeometryTypes::line.id(),double,GeometryTypes::line.id(),1>
        (result, refinementLevels(refinement), "levels");
    testStaticRefinementGeometry<GeometryTypes::line.id(),double,GeometryTypes::line.id(),1>
        (result, refinementIntervals(1<<refinement), "intervals");
  }

  // test triangle
  gt1 = gt2 = GeometryTypes::triangle;
  for (unsigned int refinement = 0; refinement < 3; refinement++)
  {
    testVirtualRefinement<double,2>(result, gt1, gt2, refinementLevels(refinement), "levels");
    testVirtualRefinement<double,2>(result, gt1, gt2, refinementIntervals(1<<refinement), "intervals");
    testStaticRefinementGeometry<GeometryTypes::triangle.id(),double,GeometryTypes::triangle.id(),2>
        (result, refinementLevels(refinement), "levels");
    testStaticRefinementGeometry<GeometryTypes::triangle.id(),double,GeometryTypes::triangle.id(),2>
        (result, refinementIntervals(1<<refinement), "intervals");
  }

  // test quadrilateral
  gt1 = gt2 = GeometryTypes::quadrilateral;
  for (unsigned int refinement = 0; refinement < 3; refinement++)
  {
    testVirtualRefinement<double,2>(result, gt1, gt2, refinementLevels(refinement), "levels");
    testVirtualRefinement<double,2>(result, gt1, gt2, refinementIntervals(1<<refinement), "intervals");
    testStaticRefinementGeometry<GeometryTypes::quadrilateral.id(),double,GeometryTypes::quadrilateral.id(),2>
        (result, refinementLevels(refinement), "levels");
    testStaticRefinementGeometry<GeometryTypes::quadrilateral.id(),double,GeometryTypes::quadrilateral.id(),2>
        (result, refinementIntervals(1<<refinement), "intervals");
  }

  // test refinement of a quadrilateral by triangles
  gt2 = GeometryTypes::triangle;
  for (unsigned int refinement = 0; refinement < 3; refinement++)
  {
    testVirtualRefinement<double,2>(result, gt1, gt2, refinementLevels(refinement), "levels");
    testVirtualRefinement<double,2>(result, gt1, gt2, refinementIntervals(1<<refinement), "intervals");
    testStaticRefinementGeometry<GeometryTypes::quadrilateral.id(),double,GeometryTypes::triangle.id(),2>
        (result, refinementLevels(refinement), "levels");
    testStaticRefinementGeometry<GeometryTypes::quadrilateral.id(),double,GeometryTypes::triangle.id(),2>
        (result, refinementIntervals(1<<refinement), "intervals");
  }

  // test tetrahedron
  gt1 = gt2 = GeometryTypes::tetrahedron;
  for (unsigned int refinement = 0; refinement < 3; refinement++)
  {
    testVirtualRefinement<double,3>(result, gt1, gt2, refinementLevels(refinement), "levels");
    testVirtualRefinement<double,3>(result, gt1, gt2, refinementIntervals(1<<refinement), "intervals");
    testStaticRefinementGeometry<GeometryTypes::tetrahedron.id(),double,GeometryTypes::tetrahedron.id(),3>
        (result, refinementLevels(refinement), "levels");
    testStaticRefinementGeometry<GeometryTypes::tetrahedron.id(),double,GeometryTypes::tetrahedron.id(),3>
        (result, refinementIntervals(1<<refinement), "intervals");
  }

  // test pyramid
  gt1 = GeometryTypes::pyramid;
  gt2 = GeometryTypes::tetrahedron;
  for (unsigned int refinement = 0; refinement < 3; refinement++)
  {
    testVirtualRefinement<double,3>(result, gt1, gt2, refinementLevels(refinement), "levels");
    testVirtualRefinement<double,3>(result, gt1, gt2, refinementIntervals(1<<refinement), "intervals");
    testStaticRefinementGeometry<GeometryTypes::pyramid.id(),double,GeometryTypes::tetrahedron.id(),3>
        (result, refinementLevels(refinement), "levels");
    testStaticRefinementGeometry<GeometryTypes::pyramid.id(),double,GeometryTypes::tetrahedron.id(),3>
        (result, refinementIntervals(1<<refinement), "intervals");
  }

  // test prism
  gt1 = GeometryTypes::prism;
  gt2 = GeometryTypes::tetrahedron;
  for (unsigned int refinement = 0; refinement < 3; refinement++)
  {
    testVirtualRefinement<double,3>(result, gt1, gt2, refinementLevels(refinement), "levels");
    testVirtualRefinement<double,3>(result, gt1, gt2, refinementIntervals(1<<refinement), "intervals");
    testStaticRefinementGeometry<GeometryTypes::prism.id(),double,GeometryTypes::tetrahedron.id(),3>
        (result, refinementLevels(refinement), "levels");
    testStaticRefinementGeometry<GeometryTypes::prism.id(),double,GeometryTypes::tetrahedron.id(),3>
        (result, refinementIntervals(1<<refinement), "intervals");
  }

  // test hexahedron
  gt1 = gt2 = GeometryTypes::hexahedron;
  for (unsigned int refinement = 0; refinement < 3; refinement++)
  {
    testVirtualRefinement<double,3>(result, gt1, gt2, refinementLevels(refinement), "levels");
    testVirtualRefinement<double,3>(result, gt1, gt2, refinementIntervals(1<<refinement), "intervals");
    testStaticRefinementGeometry<GeometryTypes::hexahedron.id(),double,GeometryTypes::hexahedron.id(),3>
        (result, refinementLevels(refinement), "levels");
    testStaticRefinementGeometry<GeometryTypes::hexahedron.id(),double,GeometryTypes::hexahedron.id(),3>
        (result, refinementIntervals(1<<refinement), "intervals");
  }

  // test refinement of hexahedron by tetrahedra
  gt1 = gt2 = GeometryTypes::hexahedron;
  for (unsigned int refinement = 0; refinement < 3; refinement++)
  {
    testVirtualRefinement<double,3>(result, gt1, gt2, refinementLevels(refinement), "levels");
    testVirtualRefinement<double,3>(result, gt1, gt2, refinementIntervals(1<<refinement), "intervals");
    testStaticRefinementGeometry<GeometryTypes::hexahedron.id(),double,GeometryTypes::tetrahedron.id(),3>
        (result, refinementLevels(refinement), "levels");
    testStaticRefinementGeometry<GeometryTypes::hexahedron.id(),double,GeometryTypes::tetrahedron.id(),3>
        (result, refinementIntervals(1<<refinement), "intervals");
  }

  return result;

}
catch (Exception &e)
{
  std::cerr << e << std::endl;
  throw;
} catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  throw;
}
