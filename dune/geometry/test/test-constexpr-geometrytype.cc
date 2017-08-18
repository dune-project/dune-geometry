#include <config.h>

#include <iostream>

#include <dune/geometry/type.hh>

int main ( int argc, char **argv )
{
  constexpr auto gt_none_1 = Dune::GeometryType();
  constexpr auto gt_none_2 = Dune::GeometryType::None();
  return not std::integral_constant<bool,gt_none_1.isNone() and gt_none_1 == gt_none_2>{};
}
