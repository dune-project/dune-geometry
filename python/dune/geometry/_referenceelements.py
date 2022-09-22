# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt
def module(dim):
    typeName = "Dune::Geo::ReferenceElement<Dune::Geo::ReferenceElementImplementation<double," + str(dim) + "> >"
    includes = ["dune/python/geometry/referenceelements.hh"]
    typeHash = "referenceelements_" + hashIt(typeName)
    generator = SimpleGenerator("ReferenceElements", "Dune::Python")
    m = generator.load(includes, typeName, typeHash)
    return m

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
