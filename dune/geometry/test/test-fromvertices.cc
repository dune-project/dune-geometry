// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <string>
#include <iostream>
#include <vector>

#include <dune/geometry/type.hh>

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
  // if (v > vertices-dim+d)
  //     return;
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
  if (vertices < dim+1 || vertices > 1<<dim)
    DUNE_THROW(Dune::Exception, "IMPOSSIBLE");
  guessTopologyId(dim,vertices,2,1,1,ids);
  if (ids.size() == 0)
    DUNE_THROW(Dune::Exception, "Impossible setting");
  std::cout << "(";
  for (size_t i=0; i<ids.size(); i++)
    std::cout << " possibility " << i << ": " << convBase(ids[i],2);
  std::cout << ") ";
  if (ids.size() > 1)
    DUNE_THROW(Dune::Exception, "Too many options");
  return ids[0];
}

void testGuess(unsigned int dim, unsigned int vertices)
{
  std::cout << "dim: " << dim;
  std::cout << " vertices: " << vertices << " ";
  unsigned int id = guessTopologyId(dim, vertices);
  std::cout << "guess: " << convBase(id, 2);
  if (dim <= 3) {
    Dune::GeometryType gt;
    gt.makeFromVertices(dim, vertices);
    std::cout << " real:  " << convBase(gt.id(), 2);
  }
  std::cout << std::endl;
}

int main()
{
  for (int d=0; d<=8; d++)
    for (int v=d+1; v<=(1<<d); v++)
    {
      try {
        testGuess(d,v);
      }
      catch (Dune::Exception & e)
      {
        std::cout << "Error: " << e.what() << std::endl;
      }
    }
}
