// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#include <config.h>

#include <limits>
#include <cmath>
#include <type_traits>
#include <vector>

#include <dune/geometry/localfiniteelementgeometry.hh>
#include <dune/geometry/mappedgeometry.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/test/checkgeometry.hh>
#include <dune/geometry/test/comparegeometries.hh>
#include <dune/geometry/test/localfiniteelements.hh>
#include <dune/geometry/test/referenceelementgeometry.hh>

template <class ctype, int cdim, Dune::GeometryType::Id id>
bool checkLocalFiniteElementGeometry ()
{
  bool pass = true;

  constexpr auto gt = Dune::GeometryType{id};
  auto refElem = Dune::referenceElement<ctype,gt.dim()>(gt);

  using LFE = std::conditional_t<
    gt.isSimplex(), Dune::Impl::P1LocalFiniteElement<ctype,ctype,gt.dim()>, std::conditional_t<
    gt.isCube(),    Dune::Impl::Q1LocalFiniteElement<ctype,ctype,gt.dim()>, void>
    >;
  auto lfe = LFE{};

  Dune::FieldMatrix<ctype,cdim,gt.dim()> A;
  for (int i = 0; i < cdim; ++i)
    for (int j = 0; j < int(gt.dim()); ++j)
      A[i][j] = ctype(3-std::abs(i-j));

  Dune::FieldVector<ctype,cdim> b;
  for (std::size_t i = 0; i < cdim; ++i)
    b[i] = ctype(i+1);

  // mapping to generate coordinates from reference-element corners
  auto f = [&](Dune::FieldVector<ctype,gt.dim()> const& x) {
    Dune::FieldVector<ctype,cdim> y;
    A.mv(x,y);
    y += b;
    return y;
  };

  auto corners = std::vector<Dune::FieldVector<ctype,cdim>>(refElem.size(gt.dim()));
  for (int i = 0; i < refElem.size(gt.dim()); ++i)
    corners[i] = f(refElem.position(i, gt.dim()));

  using Geometry = Dune::LocalFiniteElementGeometry<LFE,cdim>;

  // construct a geometry using given parametrization coefficients
  auto geometry = Geometry{refElem, lfe, corners};
  pass &= checkGeometry(geometry);

  // construct a geometry using local interpolation
  auto geometry2 = Geometry{refElem, lfe, f};
  pass &= checkGeometry(geometry2);

  // compare against MultiLinearGeometry
  using MLGeometry = Dune::MultiLinearGeometry<ctype,gt.dim(),cdim>;
  auto mlgeometry = MLGeometry{refElem, corners};
  pass &= Dune::compareGeometries(geometry, mlgeometry);

  // compare against a MappedGeometry
  using Mapping = Dune::Impl::LocalFiniteElementFunction<LFE,cdim,ctype>;
  using RefGeo = Dune::Impl::ReferenceElementGeometry<decltype(refElem)>;
  using MappedGeometry = Dune::MappedGeometry<Mapping,RefGeo>;
  auto mapping = Mapping{lfe,corners};
  auto refGeo = RefGeo{refElem};
  auto mappedgeometry = MappedGeometry{mapping, refGeo, true};
  pass &= Dune::compareGeometries(geometry, mappedgeometry);

  return pass;
}


template <class ctype>
static bool checkLocalFiniteElementGeometry ()
{
  bool pass = true;

  // pass &= checkLocalFiniteElementGeometry<ctype, 0, 0, Dune::GeometryTypes::simplex(0)>();

  pass &= checkLocalFiniteElementGeometry<ctype, 1, Dune::GeometryTypes::simplex(1)>();
  pass &= checkLocalFiniteElementGeometry<ctype, 2, Dune::GeometryTypes::simplex(1)>();
  pass &= checkLocalFiniteElementGeometry<ctype, 3, Dune::GeometryTypes::simplex(1)>();
  pass &= checkLocalFiniteElementGeometry<ctype, 4, Dune::GeometryTypes::simplex(1)>();

  pass &= checkLocalFiniteElementGeometry<ctype, 1, Dune::GeometryTypes::cube(1)>();
  pass &= checkLocalFiniteElementGeometry<ctype, 2, Dune::GeometryTypes::cube(1)>();
  pass &= checkLocalFiniteElementGeometry<ctype, 3, Dune::GeometryTypes::cube(1)>();
  pass &= checkLocalFiniteElementGeometry<ctype, 4, Dune::GeometryTypes::cube(1)>();

  pass &= checkLocalFiniteElementGeometry<ctype, 2, Dune::GeometryTypes::simplex(2)>();
  pass &= checkLocalFiniteElementGeometry<ctype, 3, Dune::GeometryTypes::simplex(2)>();
  pass &= checkLocalFiniteElementGeometry<ctype, 4, Dune::GeometryTypes::simplex(2)>();

  pass &= checkLocalFiniteElementGeometry<ctype, 2, Dune::GeometryTypes::cube(2)>();
  pass &= checkLocalFiniteElementGeometry<ctype, 3, Dune::GeometryTypes::cube(2)>();
  pass &= checkLocalFiniteElementGeometry<ctype, 4, Dune::GeometryTypes::cube(2)>();

  pass &= checkLocalFiniteElementGeometry<ctype, 3, Dune::GeometryTypes::simplex(3)>();
  pass &= checkLocalFiniteElementGeometry<ctype, 4, Dune::GeometryTypes::simplex(3)>();

  pass &= checkLocalFiniteElementGeometry<ctype, 3, Dune::GeometryTypes::cube(3)>();
  pass &= checkLocalFiniteElementGeometry<ctype, 4, Dune::GeometryTypes::cube(3)>();

  // pass &= checkLocalFiniteElementGeometry<ctype, 3, 3, Dune::GeometryTypes::pyramid>();
  // pass &= checkLocalFiniteElementGeometry<ctype, 3, 4, Dune::GeometryTypes::pyramid>();

  // pass &= checkLocalFiniteElementGeometry<ctype, 3, 3, Dune::GeometryTypes::prism>();
  // pass &= checkLocalFiniteElementGeometry<ctype, 3, 4, Dune::GeometryTypes::prism>();

  pass &= checkLocalFiniteElementGeometry<ctype, 4, Dune::GeometryTypes::simplex(4)>();
  pass &= checkLocalFiniteElementGeometry<ctype, 5, Dune::GeometryTypes::simplex(4)>();

  pass &= checkLocalFiniteElementGeometry<ctype, 4, Dune::GeometryTypes::cube(4)>();
  pass &= checkLocalFiniteElementGeometry<ctype, 5, Dune::GeometryTypes::cube(4)>();

  return pass;
}

int main ( int argc, char **argv )
{
  bool pass = true;

  std::cout << ">>> Checking ctype = double" << std::endl;
  pass &= checkLocalFiniteElementGeometry< double >();
  std::cout << ">>> Checking ctype = float" << std::endl;
  pass &= checkLocalFiniteElementGeometry< float >();

  if (!pass)
    std::cerr << "test failed!" << std::endl;

  return (pass ? 0 : 1);
}
