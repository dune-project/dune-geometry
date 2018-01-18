#include <config.h>

#include <dune/geometry/type.hh>

template<Dune::GeometryType::Id gtid>
struct Foo
{
  static constexpr Dune::GeometryType gt = gtid;
};

int main(int argc, char** argv)
{

  // make sure we can correctly roundtrip between GeometryType
  // and its Id in constexpr context
  constexpr Dune::GeometryType gt1 = Dune::GeometryTypes::triangle;

  Foo<gt1> foo;

  constexpr Dune::GeometryType gt2 = foo.gt;

  static_assert(gt1 == gt2, "The two geometry types have to compare equal");

  return 0;

}
