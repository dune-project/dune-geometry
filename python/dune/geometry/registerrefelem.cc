// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#define INCLUDE_REFELEM_INLINE
#include "registerrefelem.hh"
#ifndef DIM
#error "DIM needs to be defined!"
#endif
template void registerReferenceElementSubModule<DIM>(pybind11::module);
