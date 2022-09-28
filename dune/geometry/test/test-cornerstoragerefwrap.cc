// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

/** @file
 *
 * \brief Test std::reference_wrapper<SomeContainer> as the CornerStorage in
 *        MultiLinearGeometry.
 */

#include <config.h>

#include <iostream>
#include <functional>
#include <limits>
#include <vector>

#include <dune/common/fvector.hh>

#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/referenceelements.hh>

template<class ct>
struct TestGeometryTraits :
  public Dune::MultiLinearGeometryTraits<ct>
{
  template< int mydim, int cdim >
  struct CornerStorage
  {
    typedef std::reference_wrapper<
      const std::vector< Dune::FieldVector< ct, cdim > > > Type;
  };
};

template<class ct, int mydim, int cdim>
using TestGeometry =
  Dune::MultiLinearGeometry<ct, mydim, cdim, TestGeometryTraits<ct> >;

//! make a geometry that covers a kite shape
/**
 * \return The returned geometry will reference a static array that hold the
 *         corner information.  Because the order of destruction of static
 *         variables is a bit unclear, the returned value becomes invalid when
 *         main ends and shall not be used anymore.
 */
template<class ctype>
TestGeometry<ctype, 2, 2> kiteGeometry()
{
  static const std::vector<Dune::FieldVector<ctype, 2> > kiteCorners{
    {  0, -1.5 }, // bottom
    {  1,  0   }, // right
    { -1,  0   }, // left
    {  0,  0.5 }, // top
  };

  return { Dune::ReferenceElements<ctype, 2>::cube(), kiteCorners };
}

template<class ctype, int dim>
bool expectCenter(const Dune::FieldVector<ctype, dim> &actual,
                  const Dune::FieldVector<ctype, dim> &expected)
{
  bool match =
    (actual - expected).two_norm2() < std::numeric_limits<ctype>::epsilon();
  std::cout << (match ? "pass: " : "fail: ")
            << " Got: (" << actual << "), expected: (" << expected << ")"
            << std::endl;
  return match;
}

int main()
{
  bool pass = true;

  std::cout << "making a geometry of a kite..." << std::endl;
  auto geo = kiteGeometry<double>();
  std::cout << "checking center of kite..." << std::endl;
  pass &= expectCenter(geo.center(), { 0, -0.25 });

  {
    // the corners of an upward-pointing triangle
    std::vector<Dune::FieldVector<double, 2> > triangleCorners{
      { -1, 0 }, // left
      {  1, 0 }, // right
      {  0, 1 }, // top
    };

    std::cout << "turn the geometry into an upward-pointing triangle..."
              << std::endl;
    geo = {
      Dune::ReferenceElements<double, 2>::simplex(),
      triangleCorners,
    };
    std::cout << "checking center of upward-pointing triangle..." << std::endl;
    pass &= expectCenter(geo.center(), { 0, 1.0/3 });

    std::cout << "turning the geometry into a a leftward-pointing triangle "
              << "by moving the right corner to the bottom..." << std::endl
              << "(this is to show that the geometry really has a reference "
              << "to the coordinates, not a copy of them)" << std::endl;
    triangleCorners[1] = { 0, -1 }; // move right corner to bottom

    // NOTE: strictly speaking, the above modification invalidates geo.
    // However, we're using MultiLinearGeometry here, which does cope with
    // movement of vertices.  We don't give any guarantee that it will
    // continue to cope in the future, though.  CachedMultiLinearGeometry
    // already does not always handle this case, because it has no way of
    // knowing that it would need to update its caches.
    std::cout << "checking center of leftward-pointing triangle..."
              << std::endl;
    pass &= expectCenter(geo.center(), { -1.0/3, 0 });
    std::cout << "invalidating the geometry by letting the currently "
              << "referenced coordinate storage go out of scope..."
              << std::endl;
  }

  // triangleCorner is out of scope, so geo references freed storage and is
  // now really invalid.  The only remaining operations are assignment and
  // destruction.

  std::cout << "revalidating geometry by again assigning the geometry of a "
            << "kite..." << std::endl;
  geo = kiteGeometry<double>();
  std::cout << "checking center of kite (again)..." << std::endl;
  pass &= expectCenter(geo.center(), { 0, -0.25 });

  return pass ? 0 : 1;
}
