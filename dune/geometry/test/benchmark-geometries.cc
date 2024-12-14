// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#include <config.h>

#include <limits>
#include <cmath>
#include <type_traits>
#include <vector>

#include <dune/common/timer.hh>
#include <dune/geometry/localfiniteelementgeometry.hh>
#include <dune/geometry/mappedgeometry.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/test/checkgeometry.hh>
#include <dune/geometry/test/comparegeometries.hh>
#include <dune/geometry/test/localfiniteelements.hh>
#include <dune/geometry/test/referenceelementgeometry.hh>


template <class Geometry>
bool benchmarkGeometry (const Geometry& geo)
{
  bool pass = true;
  bool isAffine = geo.affine();
  Dune::GeometryType t = geo.type();

  const auto& quadrature = Dune::QuadratureRules<typename Geometry::ctype, Geometry::mydimension>::rule(geo.type(), 8);
  for (auto&& [pos,weight] : quadrature)
  {
    pass &= (geo.affine() == isAffine);
    pass &= (geo.type() == t);
    pass &= (geo.corners() > 0);

    for (int i = 0; i < geo.corners(); ++i) {
      pass &= ((geo.corner(i) - geo.corner((i+1)%geo.corners())).two_norm() > 0);
      pass &= ((geo.corner(i) - geo.center()).two_norm() > 0);
    }

    pass &= (geo.volume() > 0);
    pass &= (geo.global(pos).size() == Geometry::coorddimension);
    pass &= (geo.jacobian(pos).M() == Geometry::mydimension);
    pass &= (geo.jacobianTransposed(pos).N() == Geometry::mydimension);
    pass &= (geo.jacobianInverse(pos).M() == Geometry::coorddimension);
    pass &= (geo.jacobianInverseTransposed(pos).N() == Geometry::coorddimension);
    pass &= (geo.integrationElement(pos) > 0);
  }

  return pass;
}


template <class ctype, int cdim, Dune::GeometryType::Id id>
bool benchmarkGeometries (int nIter = 100)
{
  bool pass = true;

  constexpr auto gt = Dune::GeometryType{id};
  auto refElem = Dune::referenceElement<ctype,gt.dim()>(gt);
  std::cout << "Time(" << refElem.type() << "):" << std::endl;

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

  Dune::Timer t;

  // construct a geometry using given parametrization coefficients
  using Geometry = Dune::LocalFiniteElementGeometry<LFE,cdim>;
  auto geometry = Geometry{refElem, lfe, corners};
  t.reset();
  for (int i = 0; i < nIter; ++i)
    pass &= benchmarkGeometry(geometry);
  std::cout << "  LocalFiniteElementGeometry = " << t.elapsed() << "sec" << std::endl;

  // compare against MultiLinearGeometry
  using MLGeometry = Dune::MultiLinearGeometry<ctype,gt.dim(),cdim>;
  auto mlgeometry = MLGeometry{refElem, corners};
  t.reset();
  for (int i = 0; i < nIter; ++i)
    pass &= benchmarkGeometry(mlgeometry);
  std::cout << "  MultiLinearGeometry = " << t.elapsed() << "sec" << std::endl;

  // compare against a MappedGeometry
  using Mapping = Dune::Impl::LocalFiniteElementFunction<LFE,cdim,ctype>;
  using RefGeo = Dune::Impl::ReferenceElementGeometry<decltype(refElem)>;
  using MappedGeometry = Dune::MappedGeometry<Mapping,RefGeo>;
  auto mapping = Mapping{lfe,corners};
  auto refGeo = RefGeo{refElem};
  bool affine = lfe.localBasis().order() == 1 && refElem.template geometry<0>(0).affine();
  auto mappedgeometry = MappedGeometry{mapping, refGeo, affine};
  t.reset();
  for (int i = 0; i < nIter; ++i)
    pass &= benchmarkGeometry(mappedgeometry);
  std::cout << "  MappedGeometry = " << t.elapsed() << "sec" << std::endl;

  return pass;
}


template <class ctype>
static bool benchmarkGeometries (int nIter)
{
  bool pass = true;

  // pass &= benchmarkGeometries<ctype, 0, 0, Dune::GeometryTypes::simplex(0)>();

  pass &= benchmarkGeometries<ctype, 1, Dune::GeometryTypes::simplex(1)>(nIter);
  pass &= benchmarkGeometries<ctype, 2, Dune::GeometryTypes::simplex(1)>(nIter);
  pass &= benchmarkGeometries<ctype, 3, Dune::GeometryTypes::simplex(1)>(nIter);
  pass &= benchmarkGeometries<ctype, 4, Dune::GeometryTypes::simplex(1)>(nIter);

  pass &= benchmarkGeometries<ctype, 1, Dune::GeometryTypes::cube(1)>(nIter);
  pass &= benchmarkGeometries<ctype, 2, Dune::GeometryTypes::cube(1)>(nIter);
  pass &= benchmarkGeometries<ctype, 3, Dune::GeometryTypes::cube(1)>(nIter);
  pass &= benchmarkGeometries<ctype, 4, Dune::GeometryTypes::cube(1)>(nIter);

  pass &= benchmarkGeometries<ctype, 2, Dune::GeometryTypes::simplex(2)>(nIter);
  pass &= benchmarkGeometries<ctype, 3, Dune::GeometryTypes::simplex(2)>(nIter);
  pass &= benchmarkGeometries<ctype, 4, Dune::GeometryTypes::simplex(2)>(nIter);

  pass &= benchmarkGeometries<ctype, 2, Dune::GeometryTypes::cube(2)>(nIter);
  pass &= benchmarkGeometries<ctype, 3, Dune::GeometryTypes::cube(2)>(nIter);
  pass &= benchmarkGeometries<ctype, 4, Dune::GeometryTypes::cube(2)>(nIter);

  pass &= benchmarkGeometries<ctype, 3, Dune::GeometryTypes::simplex(3)>(nIter);
  pass &= benchmarkGeometries<ctype, 4, Dune::GeometryTypes::simplex(3)>(nIter);

  pass &= benchmarkGeometries<ctype, 3, Dune::GeometryTypes::cube(3)>(nIter);
  pass &= benchmarkGeometries<ctype, 4, Dune::GeometryTypes::cube(3)>(nIter);

  // pass &= benchmarkGeometries<ctype, 3, 3, Dune::GeometryTypes::pyramid>(nIter);
  // pass &= benchmarkGeometries<ctype, 3, 4, Dune::GeometryTypes::pyramid>(nIter);

  // pass &= benchmarkGeometries<ctype, 3, 3, Dune::GeometryTypes::prism>(nIter);
  // pass &= benchmarkGeometries<ctype, 3, 4, Dune::GeometryTypes::prism>(nIter);

  pass &= benchmarkGeometries<ctype, 4, Dune::GeometryTypes::simplex(4)>(nIter);
  pass &= benchmarkGeometries<ctype, 5, Dune::GeometryTypes::simplex(4)>(nIter);

  pass &= benchmarkGeometries<ctype, 4, Dune::GeometryTypes::cube(4)>(nIter);
  pass &= benchmarkGeometries<ctype, 5, Dune::GeometryTypes::cube(4)>(nIter);

  return pass;
}

int main ( int argc, char **argv )
{
  bool pass = true;

  int nIter = 10;
  if (argc > 1)
    nIter = std::atoi(argv[1]);

  std::cout << ">>> Checking ctype = double" << std::endl;
  pass &= benchmarkGeometries< double >(nIter);
  std::cout << ">>> Checking ctype = float" << std::endl;
  pass &= benchmarkGeometries< float >(nIter);

  if (!pass)
    std::cerr << "test failed!" << std::endl;

  return (pass ? 0 : 1);
}
