# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

import logging

logger = logging.getLogger(__name__)

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
    _cache = {}
    class QPQuadPyQuadrature:
        def __init__(self,quad,order,method,vertices,transform):
            self.quad_ = quad
            self.order_ = order
            self.method_ = method
            self.vertices_ = numpy.array(vertices)
            # here an error will occur if the method is invalid - how to catch?
            self.quadrature_ = getattr(quad,method)
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
        def get(self):
            return self.points_, self.weights_
        def apply(self,entity,f):
            ie = entity.geometry.integrationElement
            f_e = lambda x: f(entity,x)*ie(x)
            return self.quad_.integrate(f_e,self.vertices_,self.quadrature_)
        def __iter__(self):
            return self.quadPoints_.__iter__()
        @property
        def order(self):
            return self.quadrature_.degree
        def name(self):
            return self.method_
    def rule(gt,quadDescription):
        try:
            gt = gt.type
        except AttributeError:
            pass
        try:
            quad = _cache[(gt,quadDescription)]
            return quad
        except KeyError:
            pass
        order, method = quadDescription
        dim = gt.dim
        if gt.isLine:
            vertices = [0,1]
            quad = qp.line_segment
            from quadpy.ncube import transform
        elif gt.isTriangle:
            vertices = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]
            quad = qp.triangle
            from quadpy.nsimplex import transform
        elif gt.isQuadrilateral:
            vertices = qp.quadrilateral.rectangle_points([0.0, 1.0], [0.0, 1.0])
            quad = qp.quadrilateral
            from quadpy.ncube import transform
        elif gt.isTetrahedron:
            vertices = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
            quad = qp.tetrahedron
            from quadpy.nsimplex import transform
        elif gt.isHexahedron:
            vertices = qp.hexahedron.cube_points([0.0, 1.0], [0.0, 1.0], [0.0, 1.0])
            quad = qp.hexahedron
            from quadpy.ncube import transform
        elif gt.isPyramid:
            vertices = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
            quad = qp.pyramid
            raise ValueError("prism quadratures not yet fully supported")
        elif gt.isPrism:
            vertices = [
                [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0]
                       ]
            quad = qp.wedge
            raise ValueError("prism quadratures not yet fully supported")
        else:
            raise ValueError("no quadpy quadrature available for the geometryType " + str(gt))
        if not method and not order:
            return quad
        ret = QPQuadPyQuadrature(quad,order,method,vertices,transform)
        _cache[(gt,quadDescription)] = ret
        return ret

    def rules(methods):
        def r(entity):
            try:
                return rule(entity.type, methods[entity.type])
            except AttributeError:
                return rule(entity, methods[entity])
        return r

except ImportError as e:
    logger.warning('Unable to import quadpy: ' + " ".join(str(e).splitlines()))
    raise ImportError("Unable to import quadpy module")
