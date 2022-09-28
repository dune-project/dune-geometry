# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

import time, math, numpy
import dune.geometry as geo

def monomial(p):
    def function(point):
        return sum( x**p for x in point)
    return function

result = {3:   # integral for sum_i x_i^p over reference element
           {geo.line: 1./4.,
            geo.triangle: 0.1,
            geo.quadrilateral: 1./2.,
            geo.tetrahedron: 1./40.,
            geo.pyramid: None,
            geo.prism: None,
            geo.hexahedron: 3./4.,
           },
          4:
           {geo.line: 1./5.,
            geo.triangle: 1./15.,
            geo.quadrilateral: 2./5.,
            geo.tetrahedron: 1./70.,
            geo.pyramid: None,
            geo.prism: None,
            geo.hexahedron: 3./5.,
           },
         }



for order in [3,4]:
    rules = geo.quadratureRules(order)
    p = monomial(order)
    for t in (geo.line, geo.triangle, geo.quadrilateral,
              geo.tetrahedron, geo.pyramid, geo.prism, geo.hexahedron):
        value1 = 0
        for q in rules(t):
            value1 += p(q.position)*q.weight
        hatxs, hatws = rules(t).get()
        value2 = numpy.sum(p(hatxs) * hatws, axis=-1)
        # print(order,t,value2)
        assert abs(value1-value2)<1e-14
        if result[order][t] is not None:
            assert abs(result[order][t] - value1)<1e-12
