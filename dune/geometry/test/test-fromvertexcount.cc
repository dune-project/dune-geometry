// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <string>
#include <iostream>
#include <vector>

#include <dune/geometry/type.hh>

#include <dune/geometry/utility/typefromvertexcount.hh>

std::string convBase(unsigned long v, long base)
{
  const char* digits = "0123456789abcdef";
  std::string result;
  if((base < 2) || (base > 16)) {
    result = "Error: base out of range.";
  }
  else {
    do {
      result = digits[v % base] + result;
      v /= base;
    }
    while(v);
  }
  return result;
}

void guessTopologyId(unsigned int dim, unsigned int vertices,
                     unsigned int v, unsigned int id, unsigned int d,
                     std::vector<unsigned int> & ids)
{
  if (d == dim)
  {
    if (v == vertices)
      ids.push_back(id);
    return;
  }
  // try simplex
  guessTopologyId(dim,vertices,v+1,id,d+1,ids);
  // try cube
  guessTopologyId(dim,vertices,v*2,id+(1<<d),d+1,ids);
}

unsigned int guessTopologyId(unsigned int dim, unsigned int vertices)
{
  std::vector<unsigned int> ids;
  if (dim == 0 || dim == 1)
    return 0;
  if (vertices < dim+1 || vertices > 1u<<dim)
    DUNE_THROW(Dune::Exception, "IMPOSSIBLE");
  guessTopologyId(dim,vertices,2,1,1,ids);
  if (ids.size() == 0)
    DUNE_THROW(Dune::Exception, "Impossible setting");
  if (ids.size() > 1)
    DUNE_THROW(Dune::Exception, "Too many options");
  return ids[0];
}

void testGuess(unsigned int dim, unsigned int vertices)
{
  std::cout << "check dim: " << dim
            << " vertices: " << vertices
            << std::endl;
  unsigned int id = guessTopologyId(dim, vertices);
  Dune::GeometryType gt = Dune::geometryTypeFromVertexCount(dim, vertices);
  if (Dune::GeometryType(id,dim) != gt)
    DUNE_THROW(Dune::Exception, "Failed to guess the geometry type from the number of vertices.");
}

int main()
try {
  std::vector<std::vector<int>> configurations = { {1}, {2}, {3,4}, {4,5,6,8} };
  for (int d=0; d<=3; d++)
    for (int v : configurations[d])
        testGuess(d,v);
}
catch (Dune::Exception & e)
{
  std::cout << "Error: " << e.what() << std::endl;
  return 1;
}
