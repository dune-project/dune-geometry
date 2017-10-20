from ..generator.generator import SimpleGenerator
from dune.common.hashit import hashIt
def module(dim):
    typeName = "Dune::Geo::ReferenceElement<Dune::Geo::ReferenceElementImplementation<double," + str(dim) + "> >"
    includes = ["dune/python/geometry/referenceelements.hh"]
    typeHash = "referenceelements_" + hashIt(typeName)
    generator = SimpleGenerator("ReferenceElements", "Dune::Python")
    m = generator.load(includes, typeName, typeHash)
    return m

def referenceElement(geometryType):
    try:
        geometryType = geometryType.type
    except:
        pass
    return module(geometryType.dim).general(geometryType)
