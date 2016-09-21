#include <array>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/multilineargeometry.hh>


#include "staticIf.hh"

#define CALC_PRECISION  80
#define PRINT_PRECISION 35

using boost::multiprecision::cpp_dec_float;
using mp_type = boost::multiprecision::number<cpp_dec_float<CALC_PRECISION> >;

template <int DIM, Dune::GeometryType::BasicType>
struct IsAffine {
  static constexpr bool value = false;
};

template <int DIM>
struct IsAffine<DIM, Dune::GeometryType::simplex> {
  static constexpr bool value = true;
};

template <>
struct IsAffine<3, Dune::GeometryType::pyramid> { // whenever the base facet is a square
  static constexpr bool value = true;
};

template <>
struct IsAffine<1, Dune::GeometryType::cube> {
  static constexpr bool value = true;
};


/* Geometry descriptions:
 *
 * a0 a1 a2
 * b0 b1 b2
 * c0 c1 c2
 * ...
 **/
template <int DIM, Dune::GeometryType::BasicType t>
auto read_geometry(std::string filename)
{
  using T = mp_type;
  using LocalCoord = Dune::FieldVector<T, DIM>;
  std::ifstream in(filename, std::ios_base::in);

  std::string x_, y_, z_;
  std::vector<LocalCoord> Q;
  if (DIM == 1) {
    for (; in >> x_; )
      Q.push_back({T(x_)});
  } else if (DIM == 2) {
    for (; in >> x_ >> y_; )
      Q.push_back({T(x_), T(y_)});
  } else if (DIM == 3) {
    for (; in >> x_ >> y_ >> z_; )
      Q.push_back({T(x_), T(y_), T(z_)});
  }

  return Dune::staticIf< IsAffine<DIM,t>::value >
  (
    [&Q]() { return Dune::AffineGeometry<T,DIM,DIM>(Dune::GeometryType(t,DIM), Q); },      // then_
    [&Q]() { return Dune::MultiLinearGeometry<T,DIM,DIM>(Dune::GeometryType(t,DIM), Q); } // else_
  );
}

// coordinates given in barycentric coordinates
template <int DIM, Dune::GeometryType::BasicType t, class Out>
void print_barycentric(Out& out, std::string filename, int order,
                       std::vector<mp_type> const& data)
{
  auto geometry = read_geometry<DIM, t>(filename);

  using T = mp_type;
  using LocalCoord = Dune::FieldVector<T, DIM>;

  int corners = geometry.corners();
  int block = corners + 1;
  std::size_t npoints = data.size() / block;
  auto refelement = referenceElement(geometry);

  out << DIM << " " << order << "\n";
  for (std::size_t i = 0; i < npoints; i++) {
    LocalCoord p; p = T(0);
    for (std::size_t j = 0; j < corners; ++j)
      p += data[i*block + j] * refelement.position(j, DIM);

    for (std::size_t j = 0; j < DIM; ++j)
      out << p[j] << " ";

    out << (data[i*block + corners] * referenceElement(geometry).volume() / geometry.volume());
    out << "\n";
  }
}

template <int DIM, Dune::GeometryType::BasicType t, class Out>
void print(Out& out, std::string filename, int order,
           std::vector<mp_type> const& data)
{
  auto geometry = read_geometry<DIM, t>(filename);

  using T = mp_type;
  using LocalCoord = Dune::FieldVector<T, DIM>;

  int block = DIM + 1;
  std::size_t npoints = data.size() / block;

  out << DIM << " " << order << "\n";
  for (std::size_t i = 0; i < npoints; i++) {
    if (DIM == 1) {
      auto p = geometry.local( LocalCoord{ data[i*block] } );
      out << p[0] << " ";
    } else if (DIM == 2) {
      auto p = geometry.local( LocalCoord{ data[i*block], data[i*block+1] } );
      out << p[0] << " " << p[1] << " ";
    } else if (DIM == 3) {
      auto p = geometry.local( LocalCoord{ data[i*block], data[i*block+1], data[i*block+2] } );
      out << p[0] << " " << p[1] << " " << p[2] << " ";
    }
    out << (data[i*block + DIM] * referenceElement(geometry).volume() / geometry.volume());
    out << "\n";
  }
}


int main(int argc, char** argv)
{
  if (argc <= 3)
    throw std::runtime_error("Not enough arguments given. Usage: ./convert_rules filename order name");

  std::string filename(argv[1]);
  int order = std::atoi(argv[2]);
  std::string name(argv[3]);

  std::vector<mp_type> data;
  std::ifstream in(filename, std::ios_base::in);

  std::string x_;
  for (; in >> x_; )
    data.push_back(mp_type(x_));


  // set output precision
  std::cout << std::setprecision(PRINT_PRECISION) << std::scientific;

  // geometry description of reference geometry
  std::string geo_filename = "./rules/" + name + ".txt";
  if (name == "prism")
    print<3, Dune::GeometryType::prism>(std::cout, geo_filename, order, data);
  else if (name == "simplex2d")
    print<2, Dune::GeometryType::simplex>(std::cout, geo_filename, order, data);
  else if (name == "simplex3d")
    print<3, Dune::GeometryType::simplex>(std::cout, geo_filename, order, data);
  else if (name == "pyramid")
    print<3, Dune::GeometryType::pyramid>(std::cout, geo_filename, order, data);
  else if (name == "cube2d")
    print<2, Dune::GeometryType::cube>(std::cout, geo_filename, order, data);
  else if (name == "cube3d")
    print<3, Dune::GeometryType::cube>(std::cout, geo_filename, order, data);
  else
    throw std::runtime_error("Unknown name given");
}
