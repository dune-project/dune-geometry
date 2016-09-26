// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <dune/corepy/geometry/type.hh>
#include <dune/corepy/pybind11/pybind11.h>

PYBIND11_PLUGIN( _geometry )
{
  pybind11::module module( "_geometry" );

  Dune::CorePy::registerGeometryType( module );

  return module.ptr();
}
