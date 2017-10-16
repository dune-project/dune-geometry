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
    def __init__(self,rule,entity=None):
        self.rule_ = rule
        self.entity_ = entity
        self.points_ = None
        self.weights_ = None
    # should do this on the C++ side only the type of the return value of f
    # is unclear. Alternative: bring on the points/weights into python to avoid C++ callbacks
    def __call__(self,f,entity=None):
        if not entity:
            if not self.entity_:
                raise ValueError(\
                "The quadrature needs to be either"\
                " constructed with a given entity or the call method needs to"\
                " passed an entity")
            entity = self.entity_
        else:
            if not self.rule_.type == entity.type:
                raise ValueError(\
                "Quadrature constructed for a reference element of type="+self.rule_.type+\
                " but used to compute an integral on an entity of type="+entity.type)
        ie = entity.geometry.integrationElement
        if not self.points_:
            self.points_, self.weights_ = self.rule_.fill()
        return numpy.sum(f(entity(self.points_))*ie(self.points_)*self.weights_,axis=-1)
    def get():
        if not self.points_:
            self.points_, self.weights_ = self.rule_.fill()
        return self.points_, self.weights_
    def __iter__(self):
        return self.rule_.__iter__()
    @property
    def order(self):
        return self.rule_.order
def quadratureRule(geometryType, order):
    try:
        entity = geometryType
        geometryType = geometryType.type
    except AttributeError:
        entity = None
        pass
    rule = module(geometryType.dim).rule(geometryType,order)
    return DuneQuadrature(rule,entity)

try:
    import quadpy as qp
    import numpy
    class QPQuadPoint:
        def __init__(self,p,w):
            self.p_ = p
            self.w_ = w
        @property
        def position(self):
            return self.p_
        @property
        def weight(self):
            return self.w_
    class QPQuadPyQuadrature:
        def __init__(self,quad,order,method,vertices,transform,entity=None):
            self.quad_ = quad
            self.order_ = order
            self.method_ = method
            self.vertices_ = numpy.array(vertices)
            self.entity_ = entity
            # here an error will occur if the method is invalid - how to catch?
            self.quadrature_ = getattr(quad,method)(order)
            try:
                self.points_ = transform.transform(self.quadrature_.points.T, self.vertices_)
            except ValueError:
                self.points_ = transform.transform(self.quadrature_.points.T, self.vertices_.T).T
            try:
                self.weights_ = transform.get_detJ(self.quadrature_.points.T, self.vertices_)*self.quadrature_.weights
            except AttributeError:
                self.weights_ = transform.get_vol(self.vertices_)*self.quadrature_.weights
            self.quadPoints_ = [ QPQuadPoint(p,w) for p,w in zip(self.points_,self.weights_) ]
            self.points_ = self.points_.transpose().copy()
        def get():
            return self.points_, self.weights_
        def __call__(self,f,entity=None):
            if not entity:
                if not self.entity_:
                    raise ValueError(\
                    "The quadrature needs to be either"\
                    " constructed with a given entity or the call method needs to"\
                    " passed an entity")
                entity = self.entity_
            else:
                if not self.rule_.type == entity.type:
                    raise ValueError(\
                    "Quadrature constructed for a reference element of type="+self.rule_.type+\
                    " but used to compute an integral on an entity of type="+entity.type)
            ie = entity.geometry.integrationElement
            # f_e = lambda x:\
            #     [f(entity(x[:,i]))*ie(x[:,i])\
            #             for i in range(x.shape[1]) ]
            f_e = lambda x: f(entity(x))*ie(x)
            return self.quad_.integrate(f_e,self.vertices_,self.quadrature_)
        def __iter__(self):
            return self.quadPoints_.__iter__()
        @property
        def order(self):
            return self.quadrature_.degree
    def quadpy(geometryType,order=None,method=None):
        try:
            entity = geometryType
            geometryType = geometryType.type
        except AttributeError:
            entity = None
            pass
        gt = geometryType
        dim = gt.dim
        if gt.isLine:
            vertices = [0,1]
            quad = qp.line_segment
        elif gt.isTriangle:
            vertices = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]
            quad = qp.triangle
            from quadpy.nsimplex import tools as transform
        elif gt.isQuadrilateral:
            vertices = qp.quadrilateral.rectangle_points([0.0, 1.0], [0.0, 1.0])
            quad = qp.quadrilateral
            from quadpy.ncube import tools as transform
        elif gt.isTetrahedron:
            vertices = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
            quad = qp.tetrahedron
            from quadpy.nsimplex import tools as transform
        elif gt.isHexahedron:
            vertices = qp.hexahedron.cube_points([0.0, 1.0], [0.0, 1.0], [0.0, 1.0])
            quad = qp.hexahedron
            from quadpy.ncube import tools as transform
        elif gt.isPyramid:
            vertices = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0] [0.0, 0.0, 1.0]]
            quad = qp.pyramid
        elif gt.isPrism:
            vertices = [
                [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0]
                       ]
            quad = qp.wedge
        else:
            raise ValueError("no quadpy quadrature available for the geometryType " + str(gt))
        if not method and not order:
            return quad
        return QPQuadPyQuadrature(quad,order,method,vertices,transform,entity)

except ImportError:
    pass
