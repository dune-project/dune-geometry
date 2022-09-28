// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#include <config.h>

#include <iostream>

#include <dune/geometry/type.hh>

int main ( int /* argc */, char ** /* argv */ )
{
  constexpr auto gt_none_1 = Dune::GeometryType();
  constexpr auto gt_none_2 = Dune::GeometryTypes::none(0);
  return not std::integral_constant<
    bool,
    gt_none_1.isNone()
    and gt_none_1 == gt_none_2
    and Dune::GeometryType(1,1) == Dune::GeometryTypes::line
    >{};
}
