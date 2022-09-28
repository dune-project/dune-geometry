# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

from dune.geometry import *

# make sure t can be reconstructed from str(t)
for t in (vertex, line, triangle, quadrilateral, tetrahedron, pyramid, prism, hexahedron):
    assert GeometryType(str(t)) == t

for d in range(10):
    assert GeometryType(str(simplex(d))) == simplex(d)
    assert GeometryType(str(cube(d))) == cube(d)
    assert GeometryType(str(none(d))) == none(d)

# make sure simplices with special names can be constructed by general mechanism
for d, t in enumerate((vertex, line, triangle, tetrahedron)):
    assert GeometryType("simplex(" + str(d) + ")") == t
    assert GeometryType("general( 0, " + str(d) + ")") == t

# make sure cube with special names can be constructed by general mechanism
for d, t in enumerate((vertex, line, quadrilateral, hexahedron)):
    assert GeometryType("cube(" + str(d) + ")") == t
    assert GeometryType("general( " + str(2**d-1) + ", " + str(d) + ")") == t
