from ._geometry import *
from ._referenceelements import *
import numpy

def _duneIntegrate(self,entity,f):
    points,weights = self.fill()
    ie = entity.geometry.integrationElement
    return numpy.sum(f(entity(points))*\
            entity.geometry.integrationElement(points)*weights,axis=-1)

_duneQuadratureRules = {}
def quadratureRule(geometryType, order):
    try:
        rule = _duneQuadratureRules[(geometryType,order)]
    except KeyError:
        rule = module(geometryType.dim).rule(geometryType,order)
        setattr(rule.__class__,"apply",_duneIntegrate)
        _duneQuadratureRules[(geometryType,order)] = rule
    return rule
def quadratureRules(order):
    return lambda entity: quadratureRule(entity.type,order)

def integrate(rules,entity,f):
    return rules(entity).apply(entity,f)
