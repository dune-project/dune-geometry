from ..generator.generator import SimpleGenerator
from dune.common.hashit import hashIt
def module(dim):
    typeName = "Dune::Geo::ReferenceElement<Dune::Geo::ReferenceElementImplementation<double," + str(dim) + "> >"
    includes = ["dune/corepy/geometry/referenceelements.hh"]
    typeHash = "referenceelements_" + hashIt(typeName)
    generator = SimpleGenerator("ReferenceElements", "Dune::CorePy")
    m = generator.load(includes, typeName, typeHash)
    return m

def referenceElement(geometryType):
    try:
        geometryType = geometryType.type
    except:
        pass
    return module(geometryType.dim).general(geometryType)

# would be nice to have this on C++ side
class DuneQuadrature:
    def __init__(self,rule):
        self.rule_ = rule
    def __iter__(self):
        return self.rule_.__iter__()
    @property
    def order(self):
        return self.rule_.order
def quadratureRule(geometryType, order):
    try:
        geometryType = geometryType.type
    except:
        pass
    rule = module(geometryType.dim).rule(geometryType,order)
    return DuneQuadrature(rule)
