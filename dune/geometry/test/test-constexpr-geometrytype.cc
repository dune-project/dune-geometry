#include <config.h>

#include <iostream>

#include <dune/geometry/type.hh>

int main ( int argc, char **argv )
{
  constexpr auto gt_none_1 = Dune::GeometryType();
  return not std::integral_constant<bool,gt_none_1.isNone()>{};
}
