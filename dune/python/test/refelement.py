# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

from dune.geometry import referenceElement
from dune.geometry import vertex, line, triangle, quadrilateral, tetrahedron, pyramid, prism, hexahedron, none

def test(r):
    for codim in range(r.dimension+1):
        types = r.types(codim)
        if len(types) != r.size(codim):
            raise Exception("types tuple has wrong size")
        for i in range(len(types)):
            if types[i] != r.type(i, codim):
                raise Exception("types tuple has wrong content")

    for codim in range(r.dimension+1):
        positions = r.positions(codim)
        if len(positions) != r.size(codim):
            raise Exception("positions tuple has wrong size")
        for i in range(len(positions)):
            if positions[i] != r.position(i, codim):
                raise Exception("positions tuple has wrong content")

    if r.dimension > 0:
        normals = r.integrationOuterNormals
        if len(normals) != r.size(1):
            raise Exception("integrationOuterNormals has wrong size")
        for i in range(len(normals)):
            if normals[i] != r.integrationOuterNormal(i):
                raise Exception("integrationOuterNormals has wrong content")

test(referenceElement(vertex))
test(referenceElement(line))
test(referenceElement(triangle))
test(referenceElement(quadrilateral))
test(referenceElement(tetrahedron))
test(referenceElement(pyramid))
test(referenceElement(prism))
test(referenceElement(hexahedron))

for dim in range(4):
    if referenceElement(none(dim)) is not None:
        raise Exception("got reference element for geometry type none")
