// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#include <config.h>

#include <dune/geometry/type.hh>

#include <iostream>

template<Dune::GeometryType::Id gtid>
struct Foo
{
  static constexpr Dune::GeometryType gt = gtid;
  static unsigned int apply()
  {
    //return Foo<Dune::GeometryTypes::prismaticExtension(gt)>::gt.id(); //does not work for gcc < 10.2
    return Foo<Dune::GeometryTypes::prismaticExtension(gt).toId()>::gt.id();
  }
};

int main(int /* argc */, char** /* argv */)
{

  // make sure we can correctly roundtrip between GeometryType
  // and its Id in constexpr context
  constexpr Dune::GeometryType gt2a = Dune::GeometryTypes::triangle;

  Foo<gt2a> foo2;

  constexpr Dune::GeometryType gt2b = foo2.gt;

  static_assert(gt2a == gt2b, "The two geometry types have to compare equal");

  Foo<Dune::GeometryTypes::prismaticExtension(gt2b)> foo3;

  constexpr Dune::GeometryType gt3 = foo3.gt;

  static_assert(gt3 == Dune::GeometryTypes::prism, "The two geometry types have to compare equal");

  if (foo2.apply() != foo3.gt.id())
  {
    std::cerr << "The two topologyIds have to compare equal\n";
    return 1;
  }
  return 0;

}
