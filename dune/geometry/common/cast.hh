// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GEOMETRY_COMMON_CAST_HH
#define DUNE_GEOMETRY_COMMON_CAST_HH

#include <cstdlib>

#include <dune/geometry/common/quadmath.hh>

namespace Dune
{
  // forward declarations
  template< class K, int SIZE > class FieldVector;
  template< class K, int ROWS, int COLS > class FieldMatrix;

  namespace Impl
  {
    // General cast class that implements an `eval` method to cast a string to
    // a number of type `T`. Provide specializations for concrete floating-point
    // types below!
    template<typename T>
    struct Cast {
      static T eval(char const* str) { return T(str); }
    };

    template<> struct Cast<float> {
      static float eval(char const* str) { return std::strtof(str, NULL); }
    };
    template<> struct Cast<double> {
      static double eval(char const* str) { return std::strtod(str, NULL); }
    };
    template<> struct Cast<long double> {
      static long double eval(char const* str) { return std::strtold(str, NULL); }
    };

    template<typename T> struct Cast<FieldVector<T,1>> {
      static FieldVector<T,1> eval(char const* str) { return Cast<T>::eval(str); }
    };
    template<typename T> struct Cast<FieldMatrix<T,1,1>> {
      static FieldMatrix<T,1,1> eval(char const* str) { return Cast<T>::eval(str); }
    };

#if HAVE_FLOAT128
    // specialization for quadprecision floating-point type.
    template<> struct Cast<__float128> {
      static __float128 eval(char const* str) { return strtoflt128(str, NULL); }
    };
#endif
  }

  template<typename T>
  T cast(char const* str) { return Impl::Cast<T>::eval(str); }


} // end namespace Dune

#endif // DUNE_GEOMETRY_COMMON_CAST_HH
