# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

from ._geometry import *
from ._referenceelements import *
import numpy

def _duneIntegrate(self,entity,f):
    points,weights = self.get()
    try:
        ie = entity.geometry.integrationElement
    except AttributeError:
        ie = geometry.integrationElement
    return numpy.sum(f(entity,points)*ie(points)*weights,axis=-1)

_duneQuadratureRules = {}
def quadratureRule(geometryType, order):
    try:
        geometryType = geometryType.type
    except AttributeError:
        pass
    try:
        rule = _duneQuadratureRules[(geometryType,order)]
    except KeyError:
        rule = module(geometryType.dim).rule(geometryType,order)
        setattr(rule.__class__,"apply",_duneIntegrate)
        _duneQuadratureRules[(geometryType,order)] = rule
    return rule
def quadratureRules(order):
    return lambda entity: quadratureRule(entity,order)

def integrate(rules,entity,f):
    return rules(entity).apply(entity,f)
