# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt
def  _generateModule(dim):
    # if code is changed here also change code in dune/python/geometry/referenceelements.hh
    typeName = "Dune::Geo::ReferenceElement<Dune::Geo::ReferenceElementImplementation<double," + str(dim) + "> >"
    includes = ["dune/python/geometry/referenceelements.hh"]
    typeHash = "referenceelements_" + hashIt(typeName)
    generator = SimpleGenerator("ReferenceElements", "Dune::Python")
    m = generator.load(includes, typeName, typeHash)
    return m

def module(dim):
    # try to load default module added in _geometry.cc
    # this needs the DUNE_ENABLE_PYTHONMODULE_PRECOMPILE to be enabled
    # otherwise all modules will be generated and compiled
    try:
        import importlib
        return importlib.import_module( "dune.geometry._geometry._defaultreferenceelements_"+str(dim) )
    except ImportError:
        return _generateModule(dim)


_duneReferenceElements = {}
def referenceElement(geometryType):
    try:
        geometryType = geometryType.type
    except:
        pass
    try:
        ref = _duneReferenceElements[geometryType]
    except KeyError:
        ref = module(geometryType.dim).general(geometryType)
        _duneReferenceElements[geometryType] = ref
    return ref
