from ..generator.generator import SimpleGenerator
from dune.common.hashit import hashIt

def _generateModule(dim):
    # if code is changed here also change code in dune/python/geometry/referenceelements.hh
    typeName = "Dune::Geo::ReferenceElement<Dune::Geo::ReferenceElementImplementation<double," + str(dim) + "> >"
    includes = ["dune/python/geometry/referenceelements.hh"]
    typeHash = "referenceelements_" + hashIt(typeName)
    generator = SimpleGenerator("ReferenceElements", "Dune::Python")
    m = generator.load(includes, typeName, typeHash)
    return m

def module(dim):
    # try to load default module
    try:
        import importlib
        return importlib.import_module( "dune.geometry._defaultreferenceelements_"+str(dim) )
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
        # store type in container
        _duneReferenceElements[geometryType] = ref
    return ref
