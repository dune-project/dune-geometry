// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>
#define USING_DUNE_PYTHON 1
#include <dune/python/geometry/referenceelements.hh>

PYBIND11_MODULE( _defaultreferenceelements_3, module )
{
  Dune::Python::registerReferenceElementToModule< 3 >( module );
}
