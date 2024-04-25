// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#ifndef DUNE_GEOMETRY_QUADRATURERULES_NUMBERFROMSTRING_HH
#define DUNE_GEOMETRY_QUADRATURERULES_NUMBERFROMSTRING_HH

#include <type_traits>

//! expand the number to a value and to a string
#define DUNE_NUMBER_FROM_STRING(type,value) \
  Dune::Impl::numberFromString< type >(value, #value)

namespace Dune::Impl {

//! Construct the number type `ct` from `double` or from a character sequence
template<typename ct>
ct numberFromString([[maybe_unused]] double value, [[maybe_unused]] const char* str)
{
  if constexpr(std::is_constructible_v<ct,const char*>)
    return ct{str};
  else
    return value;
}

} // end namespace Dune::Impl

#endif // DUNE_GEOMETRY_QUADRATURERULES_NUMBERFROMSTRING_HH
