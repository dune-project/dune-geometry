// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

/** \file
    \brief A unit test for the AxisAlignedCubeGeometry class
 */

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/ios_state.hh>

#include <dune/geometry/axisalignedcubegeometry.hh>
#include <dune/geometry/test/checkgeometry.hh>


using namespace Dune;

void fail(int &result) {
  result = 1;
}
void pass(int &result) {
  if(result == 77) result = 0;
}

/** \brief Test whether the given BasicGeometry object 'affine' attribute is
 *         set correctly
 *
 * \param geometry       The Geometry object to test.
 * \param expectedAffine Whether the geometry should be affine.
 * \param result         Collect pass/fail results.
 */
template <class TestGeometry>
void testBasicGeometryAffine(const TestGeometry& geometry, int &result)
{
  if( not geometry.affine() ) {
    Dune::ios_base_all_saver saver(std::cerr);
    std::cerr << std::boolalpha;
    std::cerr << "Error: Hypercube geometry is not affine!" << std::endl;
    fail(result);
  }
  else
    pass(result);
}


template <int dim, int coorddim>
void testCodimZero(int& result)
{
  std::cout << "== coorddimension = " << coorddim << std::endl;
  std::cout << "=== mydimension = " << dim << std::endl;

  typedef AxisAlignedCubeGeometry<double, dim, coorddim> ElementGeometry;

  FieldVector<double,coorddim> lower(0);
  FieldVector<double,coorddim> upper(1);
  std::bitset<coorddim> axes(0);
  for (int i=0; i<dim; i++)
    axes[i] = true;

  ElementGeometry geometry( lower, upper, axes );

  if (checkGeometry(geometry))
    pass(result);
  else
    fail(result);

  testBasicGeometryAffine(geometry, result);

  // test assignability
  ElementGeometry geometry2( lower, upper, axes );
  geometry2 = geometry;
}

template <int dim, int coorddim>
void testCodimNonZero(int& result)
{
  std::cout << "== coorddimension = " << coorddim << std::endl;
  std::cout << "=== mydimension = " << dim << std::endl;

  typedef AxisAlignedCubeGeometry<double, dim, coorddim> ElementGeometry;

  FieldVector<double,coorddim> lower(-1);
  FieldVector<double,coorddim> upper(3);

  for (size_t i=0; i<(1<<coorddim)-1; i++) {

    std::bitset<coorddim> axes(i);
    if (axes.count() != dim)
      continue;

    ElementGeometry geometry( lower, upper, axes );

    if (checkGeometry(geometry))
      pass(result);
    else
      fail(result);

    testBasicGeometryAffine(geometry, result);

    // test assignability
    ElementGeometry geometry2( lower, upper, axes );
    geometry2 = geometry;
  }

}


int main (int /* argc */ , char ** /* argv */) try
{
  // 77 means "SKIP"
  int result = 77;

  testCodimZero<0,0>(result);
  testCodimZero<1,1>(result);
  testCodimZero<2,2>(result);
  testCodimZero<3,3>(result);

  testCodimNonZero<0,2>(result);
  testCodimNonZero<1,2>(result);
  testCodimNonZero<2,2>(result);

  testCodimNonZero<0,3>(result);
  testCodimNonZero<1,3>(result);
  testCodimNonZero<2,3>(result);
  testCodimNonZero<3,3>(result);

  // Test what happens when a zero-dimensional geometry (a point) is
  // constructed with the codim-zero constructor (taking lower and upper
  // corners).
  testCodimZero<0,3>(result);

  return result;
}
catch (Dune::Exception& e) {
  std::cerr << e << std::endl;
  throw;
} catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  throw;
}
