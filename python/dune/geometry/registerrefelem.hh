// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#pragma once

#include <string>
#include <dune/common/tupleutility.hh>
#include <dune/python/geometry/type.hh>
#include <dune/python/geometry/referenceelements.hh>
#include <dune/python/pybind11/pybind11.h>

template <int dim>
void registerReferenceElementSubModule(pybind11::module module)
#ifdef INCLUDE_REFELEM_INLINE
{
  // add commonly used reference elements to the library
  std::string name = std::string("_defaultreferenceelements_");
  name += std::to_string(dim);
  std::string desc = "Default Reference Elements dim=";
  desc += std::to_string(dim);
  auto refElem = module.def_submodule( name.c_str(), desc.c_str() );
  Dune::Python::registerReferenceElementToModule< dim >( refElem );
}
#else
;
#endif
