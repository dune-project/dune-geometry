// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GEOMETRY_COMMON_QUADMATH_HH
#define DUNE_GEOMETRY_COMMON_QUADMATH_HH

#if HAVE_FLOAT128
#include <ostream>
#include <quadmath.h>

namespace Dune
{
  template<typename CharT, typename Traits>
  std::basic_ostream<CahrT,Traits>& operator<< (std::basic_ostream<CahrT,Traits>& out,
                                                __float128 const& value)
  {
    const std::size_t bufSize = 128;
    CharT buf[128];

    std::string format = "%." + std::to_string(out.precision()) + "Q" +
                         (out.flags() | std::ios_base::scientific) ? "e" : "f";
    const int numChars = quadmath_snprintf(buf, bufSize, format.c_str(), value);
    if (size_t(numChars) >= bufSize) {
      DUNE_THROW(RangeError, "Failed to print __float128 value: buffer overflow");
    }
    out << buf;
    return out;
  }

} // end namespace Dune


#endif // HAVE_FLOAT128
#endif // DUNE_GEOMETRY_COMMON_QUADMATH_HH
