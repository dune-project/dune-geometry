// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
/** \file
    \brief A unit test for the ReferenceElements
    \todo For several element types the subEntities() method is not tested yet!
 */

#include <config.h>

#include <iostream>

#include <dune/geometry/referenceelements.hh>

using namespace Dune;

#define test(a) if (! (a) ) { std::cerr << __FILE__ << ":" << __LINE__ << ": Test `" # a "' failed" << std::endl; errors++; }
#define testcmp(a,b) if (! (a == b) ) { std::cerr << __FILE__ << ":" << __LINE__ << ": Test `" # a " == " # b "' failed (got " << a << ")" << std::endl; errors++; }



template<class RE>
int checkSubEntities(const RE& re)
{
  int errors = 0;

  for (std::size_t codim = 0; codim <= RE::dimension; ++codim)
  {
    for (std::size_t i = 0; i < std::size_t(re.size(codim)); ++i)
    {
      for (std::size_t c = 0; c <= RE::dimension; ++c)
      {
        auto subEntities = re.subEntities(i, codim, c);
        auto it = subEntities.begin();
        auto end = subEntities.end();

        // check size of subEntities range
        testcmp(re.size(i, codim, c), end-it);
        testcmp(std::size_t(re.size(i, codim, c)), subEntities.size());

        // check if sub-subentity range is empty if
        // sub-subentity codim c < codim of subentity
        if (c<codim)
          testcmp(re.size(i, codim, c), 0);

        // check if sub-subentity range is singleton if
        // sub-subentity codim c = codim of subentity
        if (c==codim)
          testcmp(re.size(i, codim, c), 1);

        // check if sub-subentity range has multiple
        // entries for sub-subentity codim c > codim of subentity
        if (c>codim)
          test(re.size(i, codim, c) > 1);

        // check if subEntities() is conforming with subEntity()
        for(std::size_t j = 0; j < subEntities.size(); ++j)
        {
          testcmp(*it, std::size_t(re.subEntity(i, codim, j, c)));
          ++it;
        }

        // check is contains is true for all indices in the range
        for(auto j : subEntities)
          test(subEntities.contains(j));

        // check if contains is consistent
        std::vector<bool> containedSubEntities(re.size(c), false);
        for(auto j : subEntities)
          containedSubEntities[j] = true;
        for(auto j: range(std::size_t(re.size(c))))
          testcmp(containedSubEntities[j],subEntities.contains(j));

        // after incrementing size times, the range should be exhausted
        testcmp(it, end);
      }
    }
  }
  return errors;
}



int main () try
{
  int errors = 0;

  // //////////////////////////////////////////////////////////////////////////
  //   We test the different elements separately.  I will not try to be smart
  //   and write routines that check all types at once.  That would just be
  //   illegible.
  // //////////////////////////////////////////////////////////////////////////

  GeometryType type;

  // //////////////////////////////////////////////////////////////////////////
  //   Test segment
  // //////////////////////////////////////////////////////////////////////////

  type = GeometryTypes::line;

  {
    // check default constructibility and comparison operators
    ReferenceElements<double,1>::ReferenceElement r1, r2;
    test(r1 == r2);
    test(not (r1 != r2));
    // check hash value
    testcmp(hash_value(r1),0);
  }

  const Transitional::ReferenceElement<double,Dim<1>> referenceLine = ReferenceElements<double, 1>::general(type);

  // size(int c)
  testcmp(referenceLine.size(0),1);
  testcmp(referenceLine.size(1),2);

  // size(int i, int c, int cc)
  testcmp(referenceLine.size(0,0,0),1);
  testcmp(referenceLine.size(0,0,1),2);
  testcmp(referenceLine.size(0,1,1),1);
  testcmp(referenceLine.size(1,1,1),1);

  // subEntity(int i, int c, int ii, int cc)

  // type(int i, int c)
  test(referenceLine.type(0,0).isLine());
  test(referenceLine.type(0,1).isVertex());
  test(referenceLine.type(1,1).isVertex());

  // test the 'volume' method
  decltype(referenceLine)::Volume lineVolume = referenceLine.volume();
  testcmp(lineVolume, 1);

  // test the 'geometry' method
  const ReferenceElements<double,1>::ReferenceElement::Codim<0>::Geometry referenceLineMapping = referenceLine.geometry< 0 >( 0 );
  referenceLineMapping.corner(0);

  errors += checkSubEntities(referenceLine);

  // //////////////////////////////////////////////////////////////////////////
  //   Test triangle
  // //////////////////////////////////////////////////////////////////////////

  type = GeometryTypes::triangle;

  const auto referenceTriangle = referenceElement<double,2>(type);

  // size(int c)
  testcmp(referenceTriangle.size(0),1);
  testcmp(referenceTriangle.size(1),3);
  testcmp(referenceTriangle.size(2),3);

  // size(int i, int c, int cc)
  testcmp(referenceTriangle.size(0,0,0),1);
  testcmp(referenceTriangle.size(0,0,1),3);
  testcmp(referenceTriangle.size(0,0,2),3);

  testcmp(referenceTriangle.size(0,1,1),1);
  testcmp(referenceTriangle.size(0,1,2),2);
  testcmp(referenceTriangle.size(1,1,1),1);
  testcmp(referenceTriangle.size(1,1,2),2);
  testcmp(referenceTriangle.size(2,1,1),1);
  testcmp(referenceTriangle.size(2,1,2),2);

  testcmp(referenceTriangle.size(0,2,2),1);
  testcmp(referenceTriangle.size(1,2,2),1);
  testcmp(referenceTriangle.size(2,2,2),1);

  // subEntity(int i, int c, int ii, int cc)
  testcmp(referenceTriangle.subEntity(0,0,0,0),0);
  testcmp(referenceTriangle.subEntity(0,0,0,1),0);
  testcmp(referenceTriangle.subEntity(0,0,1,1),1);
  testcmp(referenceTriangle.subEntity(0,0,2,1),2);
  testcmp(referenceTriangle.subEntity(0,0,0,2),0);
  testcmp(referenceTriangle.subEntity(0,0,1,2),1);
  testcmp(referenceTriangle.subEntity(0,0,2,2),2);

  testcmp(referenceTriangle.subEntity(0,1,0,1),0);
  testcmp(referenceTriangle.subEntity(1,1,0,1),1);
  testcmp(referenceTriangle.subEntity(2,1,0,1),2);

  testcmp(referenceTriangle.subEntity(0,1,0,2),0);
  testcmp(referenceTriangle.subEntity(0,1,1,2),1);
  testcmp(referenceTriangle.subEntity(1,1,0,2),0);
  testcmp(referenceTriangle.subEntity(1,1,1,2),2);
  testcmp(referenceTriangle.subEntity(2,1,0,2),1);
  testcmp(referenceTriangle.subEntity(2,1,1,2),2);

  testcmp(referenceTriangle.subEntity(0,2,0,2),0);
  testcmp(referenceTriangle.subEntity(1,2,0,2),1);
  testcmp(referenceTriangle.subEntity(2,2,0,2),2);

  // type(int i, int c)
  test(referenceTriangle.type(0,0).isTriangle());

  test(referenceTriangle.type(0,1).isLine());
  test(referenceTriangle.type(1,1).isLine());
  test(referenceTriangle.type(2,1).isLine());

  test(referenceTriangle.type(0,2).isVertex());
  test(referenceTriangle.type(1,2).isVertex());
  test(referenceTriangle.type(2,2).isVertex());

  // test the 'volume' method
  decltype(referenceTriangle)::Volume triangleVolume = referenceTriangle.volume();
  testcmp(triangleVolume, 0.5);

  // test the 'geometry' method
  const Transitional::ReferenceElement<double,Dim<2>>::Codim<0>::Geometry referenceTriangleMapping = referenceTriangle.geometry< 0 >( 0 );
  referenceTriangleMapping.corner(0);

  // test the checkInside method
  test(referenceTriangle.checkInside({0.3,0.3}));

  errors += checkSubEntities(referenceTriangle);

  // //////////////////////////////////////////////////////////////////////////
  //   Test quadrilateral
  // //////////////////////////////////////////////////////////////////////////

  type = GeometryTypes::quadrilateral;

  const Transitional::ReferenceElement<double,Dim<2>> referenceQuad = referenceElement<double>(type,Dim<2>());

  // size(int c)
  testcmp(referenceQuad.size(0),1);
  testcmp(referenceQuad.size(1),4);
  testcmp(referenceQuad.size(2),4);

  // size(int i, int c, int cc)
  testcmp(referenceQuad.size(0,0,0),1);
  testcmp(referenceQuad.size(0,0,1),4);
  testcmp(referenceQuad.size(0,0,2),4);

  testcmp(referenceQuad.size(0,1,1),1);
  testcmp(referenceQuad.size(0,1,2),2);
  testcmp(referenceQuad.size(1,1,1),1);
  testcmp(referenceQuad.size(1,1,2),2);
  testcmp(referenceQuad.size(2,1,1),1);
  testcmp(referenceQuad.size(2,1,2),2);
  testcmp(referenceQuad.size(3,1,1),1);
  testcmp(referenceQuad.size(3,1,2),2);

  testcmp(referenceQuad.size(0,2,2),1);
  testcmp(referenceQuad.size(1,2,2),1);
  testcmp(referenceQuad.size(2,2,2),1);
  testcmp(referenceQuad.size(3,2,2),1);

  // subEntity(int i, int c, int ii, int cc)
  testcmp(referenceQuad.subEntity(0,0,0,0),0);
  testcmp(referenceQuad.subEntity(0,0,0,1),0);
  testcmp(referenceQuad.subEntity(0,0,1,1),1);
  testcmp(referenceQuad.subEntity(0,0,2,1),2);
  testcmp(referenceQuad.subEntity(0,0,3,1),3);
  testcmp(referenceQuad.subEntity(0,0,0,2),0);
  testcmp(referenceQuad.subEntity(0,0,1,2),1);
  testcmp(referenceQuad.subEntity(0,0,2,2),2);
  testcmp(referenceQuad.subEntity(0,0,3,2),3);

  testcmp(referenceQuad.subEntity(0,1,0,1),0);
  testcmp(referenceQuad.subEntity(1,1,0,1),1);
  testcmp(referenceQuad.subEntity(2,1,0,1),2);
  testcmp(referenceQuad.subEntity(3,1,0,1),3);

  testcmp(referenceQuad.subEntity(0,1,0,2),0);
  testcmp(referenceQuad.subEntity(0,1,1,2),2);
  testcmp(referenceQuad.subEntity(1,1,0,2),1);
  testcmp(referenceQuad.subEntity(1,1,1,2),3);
  testcmp(referenceQuad.subEntity(2,1,0,2),0);
  testcmp(referenceQuad.subEntity(2,1,1,2),1);
  testcmp(referenceQuad.subEntity(3,1,0,2),2);
  testcmp(referenceQuad.subEntity(3,1,1,2),3);

  testcmp(referenceQuad.subEntity(0,2,0,2),0);
  testcmp(referenceQuad.subEntity(1,2,0,2),1);
  testcmp(referenceQuad.subEntity(2,2,0,2),2);
  testcmp(referenceQuad.subEntity(3,2,0,2),3);

  // type(int i, int c)
  test(referenceQuad.type(0,0).isQuadrilateral());

  test(referenceQuad.type(0,1).isLine());
  test(referenceQuad.type(1,1).isLine());
  test(referenceQuad.type(2,1).isLine());
  test(referenceQuad.type(3,1).isLine());

  test(referenceQuad.type(0,2).isVertex());
  test(referenceQuad.type(1,2).isVertex());
  test(referenceQuad.type(2,2).isVertex());
  test(referenceQuad.type(3,2).isVertex());

  errors += checkSubEntities(referenceQuad);

  // test the 'volume' method
  decltype(referenceQuad)::Volume quadVolume = referenceQuad.volume();
  testcmp(quadVolume, 1);

  // test the 'geometry' method
  const Transitional::ReferenceElement<double,Dim<2>>::Codim<0>::Geometry referenceQuadMapping = referenceQuad.geometry< 0 >( 0 );
  referenceQuadMapping.corner(0);

  // //////////////////////////////////////////////////////////////////////////
  //   Test tetrahedron
  // //////////////////////////////////////////////////////////////////////////

  type = GeometryTypes::tetrahedron;

  const Transitional::ReferenceElement<double,Dim<3>> referenceTetra = referenceElement(double(),type,Dim<3>());

  // size(int c)
  testcmp(referenceTetra.size(0),1);
  testcmp(referenceTetra.size(1),4);
  testcmp(referenceTetra.size(2),6);
  testcmp(referenceTetra.size(3),4);

  // size(int i, int c, int cc)
  testcmp(referenceTetra.size(0,0,0),1);
  testcmp(referenceTetra.size(0,0,1),4);
  testcmp(referenceTetra.size(0,0,2),6);
  testcmp(referenceTetra.size(0,0,3),4);

  for (int i=0; i<referenceTetra.size(1); i++) {
    testcmp(referenceTetra.size(i,1,1),1);
    testcmp(referenceTetra.size(i,1,2),3);
    testcmp(referenceTetra.size(i,1,3),3);
  }

  for (int i=0; i<referenceTetra.size(2); i++) {
    testcmp(referenceTetra.size(i,2,2),1);
    testcmp(referenceTetra.size(i,2,3),2);
  }

  for (int i=0; i<referenceTetra.size(3); i++)
    testcmp(referenceTetra.size(i,3,3),1);

  // subEntity(int i, int c, int ii, int cc)


  // type(int i, int c)
  test(referenceTetra.type(0,0).isTetrahedron());

  for (int i=0; i<referenceTetra.size(1); i++)
    test(referenceTetra.type(i,1).isTriangle());

  for (int i=0; i<referenceTetra.size(2); i++)
    test(referenceTetra.type(i,2).isLine());

  for (int i=0; i<referenceTetra.size(3); i++)
    test(referenceTetra.type(i,3).isVertex());

  // test the 'volume' method
  decltype(referenceTetra)::Volume tetraVolume = referenceTetra.volume();
  testcmp(tetraVolume, 1.0/6.0);

  // test the 'geometry' method
  const decltype(referenceTetra)::Codim<0>::Geometry referenceTetraMapping = referenceTetra.geometry< 0 >( 0 );
  referenceTetraMapping.corner(0);

  // test the checkInside method
  test(referenceTetra.checkInside({0.3,0.3,0.3}));

  errors += checkSubEntities(referenceTetra);

  // //////////////////////////////////////////////////////////////////////////
  //   Test pyramid
  // //////////////////////////////////////////////////////////////////////////

  type = GeometryTypes::pyramid;

  const auto referencePyramid = ReferenceElements<double, 3>::general(type);

  // size(int c)
  testcmp(referencePyramid.size(0),1);
  testcmp(referencePyramid.size(1),5);
  testcmp(referencePyramid.size(2),8);
  testcmp(referencePyramid.size(3),5);

  // size(int i, int c, int cc)
  testcmp(referencePyramid.size(0,0,0),1);
  testcmp(referencePyramid.size(0,0,1),5);
  testcmp(referencePyramid.size(0,0,2),8);
  testcmp(referencePyramid.size(0,0,3),5);

  testcmp(referencePyramid.size(0,1,1),1);
  testcmp(referencePyramid.size(0,1,2),4);
  testcmp(referencePyramid.size(0,1,3),4);
  testcmp(referencePyramid.size(1,1,1),1);
  testcmp(referencePyramid.size(1,1,2),3);
  testcmp(referencePyramid.size(1,1,3),3);
  testcmp(referencePyramid.size(2,1,1),1);
  testcmp(referencePyramid.size(2,1,2),3);
  testcmp(referencePyramid.size(2,1,3),3);
  testcmp(referencePyramid.size(3,1,1),1);
  testcmp(referencePyramid.size(3,1,2),3);
  testcmp(referencePyramid.size(3,1,3),3);
  testcmp(referencePyramid.size(4,1,1),1);
  testcmp(referencePyramid.size(4,1,2),3);
  testcmp(referencePyramid.size(4,1,3),3);

  for (int i=0; i<referencePyramid.size(2); i++) {
    testcmp(referencePyramid.size(i,2,2),1);
    testcmp(referencePyramid.size(i,2,3),2);
  }

  for (int i=0; i<referencePyramid.size(3); i++)
    testcmp(referencePyramid.size(i,3,3),1);

  // subEntity(int i, int c, int ii, int cc)

  // type(int i, int c)
  test(referencePyramid.type(0,0).isPyramid());

  test(referencePyramid.type(0,1).isQuadrilateral());
  test(referencePyramid.type(1,1).isTriangle());
  test(referencePyramid.type(2,1).isTriangle());
  test(referencePyramid.type(3,1).isTriangle());
  test(referencePyramid.type(4,1).isTriangle());

  for (int i=0; i<referencePyramid.size(2); i++)
    test(referencePyramid.type(i,2).isLine());

  for (int i=0; i<referencePyramid.size(3); i++)
    test(referencePyramid.type(i,3).isVertex());

  // test the 'volume' method
  decltype(referencePyramid)::Volume pyramidVolume = referencePyramid.volume();
  testcmp(pyramidVolume, 1.0/3.0);

  // test the 'geometry' method
  const auto referencePyramidMapping = referencePyramid.geometry< 0 >( 0 );
  referencePyramidMapping.corner(0);

  errors += checkSubEntities(referencePyramid);

  // //////////////////////////////////////////////////////////////////////////
  //   Test prism
  // //////////////////////////////////////////////////////////////////////////

  type = GeometryTypes::prism;

  const auto referencePrism = referenceElement<double,3>(type);

  // size(int c)
  testcmp(referencePrism.size(0),1);
  testcmp(referencePrism.size(1),5);
  testcmp(referencePrism.size(2),9);
  testcmp(referencePrism.size(3),6);

  // size(int i, int c, int cc)
  testcmp(referencePrism.size(0,0,0),1);
  testcmp(referencePrism.size(0,0,1),5);
  testcmp(referencePrism.size(0,0,2),9);
  testcmp(referencePrism.size(0,0,3),6);

  testcmp(referencePrism.size(0,1,1),1);
  testcmp(referencePrism.size(0,1,2),4);
  testcmp(referencePrism.size(0,1,3),4);
  testcmp(referencePrism.size(1,1,1),1);
  testcmp(referencePrism.size(1,1,2),4);
  testcmp(referencePrism.size(1,1,3),4);
  testcmp(referencePrism.size(2,1,1),1);
  testcmp(referencePrism.size(2,1,2),4);
  testcmp(referencePrism.size(2,1,3),4);
  testcmp(referencePrism.size(3,1,1),1);
  testcmp(referencePrism.size(3,1,2),3);
  testcmp(referencePrism.size(3,1,3),3);
  testcmp(referencePrism.size(4,1,1),1);
  testcmp(referencePrism.size(4,1,2),3);
  testcmp(referencePrism.size(4,1,3),3);

  for (int i=0; i<referencePrism.size(2); i++) {
    testcmp(referencePrism.size(i,2,2),1);
    testcmp(referencePrism.size(i,2,3),2);
  }

  for (int i=0; i<referencePrism.size(3); i++)
    testcmp(referencePrism.size(i,3,3),1);

  // subEntity(int i, int c, int ii, int cc)

  // type(int i, int c)
  test(referencePrism.type(0,0).isPrism());

  test(referencePrism.type(0,1).isQuadrilateral());
  test(referencePrism.type(1,1).isQuadrilateral());
  test(referencePrism.type(2,1).isQuadrilateral());
  test(referencePrism.type(3,1).isTriangle());
  test(referencePrism.type(4,1).isTriangle());

  for (int i=0; i<referencePrism.size(2); i++)
    test(referencePrism.type(i,2).isLine());

  for (int i=0; i<referencePrism.size(3); i++)
    test(referencePrism.type(i,3).isVertex());

  // test the 'volume' method
  decltype(referencePrism)::Volume prismVolume = referencePrism.volume();
  testcmp(prismVolume, 0.5);

  // test the 'geometry' method
  const auto referencePrismMapping = referencePrism.geometry< 0 >( 0 );
  referencePrismMapping.corner(0);

  errors += checkSubEntities(referencePrism);

  // //////////////////////////////////////////////////////////////////////////
  //   Test hexahedron
  // //////////////////////////////////////////////////////////////////////////

  type = GeometryTypes::hexahedron;

  const auto referenceHexa = ReferenceElements<double, 3>::general(type);

  // size(int c)
  testcmp(referenceHexa.size(0),1);
  testcmp(referenceHexa.size(1),6);
  testcmp(referenceHexa.size(2),12);
  testcmp(referenceHexa.size(3),8);

  // size(int i, int c, int cc)
  testcmp(referenceHexa.size(0,0,0),1);
  testcmp(referenceHexa.size(0,0,1),6);
  testcmp(referenceHexa.size(0,0,2),12);
  testcmp(referenceHexa.size(0,0,3),8);

  for (int i=0; i<referenceHexa.size(1); i++) {
    testcmp(referenceHexa.size(i,1,1),1);
    testcmp(referenceHexa.size(i,1,2),4);
    testcmp(referenceHexa.size(i,1,3),4);
  }

  for (int i=0; i<referenceHexa.size(2); i++) {
    testcmp(referenceHexa.size(i,2,2),1);
    testcmp(referenceHexa.size(i,2,3),2);
  }

  for (int i=0; i<referenceHexa.size(3); i++)
    testcmp(referenceHexa.size(i,3,3),1);

  // subEntity(int i, int c, int ii, int cc)

  // type(int i, int c)
  test(referenceHexa.type(0,0).isHexahedron());

  for (int i=0; i<referenceHexa.size(1); i++)
    test(referenceHexa.type(i,1).isQuadrilateral());

  for (int i=0; i<referenceHexa.size(2); i++)
    test(referenceHexa.type(i,2).isLine());

  for (int i=0; i<referenceHexa.size(3); i++)
    test(referenceHexa.type(i,3).isVertex());

  // test the 'volume' method
  decltype(referenceHexa)::Volume hexaVolume = referenceHexa.volume();
  testcmp(hexaVolume, 1);

  // test the 'geometry' method
  const Transitional::ReferenceElement<double,Dim<3>>::Codim<0>::Geometry referenceHexaMapping = referenceHexa.geometry< 0 >( 0 );
  referenceHexaMapping.corner(0);

  errors += checkSubEntities(referenceHexa);

  return errors>0 ? 1 : 0;

}
catch (Exception &e) {
  std::cerr << e << std::endl;
  return 1;
} catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
