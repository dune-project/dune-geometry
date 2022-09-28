// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#include <config.h>

#include <functional>
#include <vector>

#include <dune/common/fvector.hh>

#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/geometry/test/checkgeometry.hh>


template<class ct>
struct ReferenceWrapperGeometryTraits :
  public Dune::MultiLinearGeometryTraits<ct>
{
  template< int mydim, int cdim >
  struct CornerStorage
  {
    typedef std::reference_wrapper<
      const std::vector< Dune::FieldVector< ct, cdim > > > Type;
  };
};

template< class ctype, int mydim, int cdim >
static Dune::FieldVector< ctype, cdim >
map ( const Dune::FieldMatrix< ctype, mydim, mydim > &A,
      const Dune::FieldMatrix< ctype, cdim, cdim > &B,
      const Dune::FieldVector< ctype, mydim > &x )
{
  // compute Ax
  Dune::FieldVector< ctype, mydim > Ax;
  A.mv( x, Ax );

  // embed Ax into the larger space (eAx)
  Dune::FieldVector< ctype, cdim > eAx( 0 );
  for( int j = 0; j < mydim; ++j )
    eAx[ j ] = Ax[ j ];

  // combute y = B eAx
  Dune::FieldVector< ctype, cdim > y;
  B.mv( eAx, y );
  return y;
}

template <class RefElement, class Geometry, class GlobalCoordinate>
static bool testCheckInside( const RefElement& refElement,
                             const Geometry& geometry,
                             const GlobalCoordinate& global,
                             const bool testResult )
{
  auto local = geometry.local( global );
  std::cout << "Checking geometry " << geometry.type() << std::endl;
  if( refElement.checkInside( local ) != testResult )
  {
    std::cerr << "Error: Local operation failed!" << std::endl;
    return false;
  }
  else if( testResult ) // for positive result check global and local
  {
    // double check that global and local are inverse to each other
    auto globalCheck = geometry.global( local );
    if( (globalCheck - global).two_norm() > 1e-8 )
    {
      std::cerr << " global and local not inverse to each other!" << std::endl;
      std::cerr << "(" << globalCheck << ") is not (" << global << ")" << std::endl;
      return false ;
    }
  }
  std::string name = testResult ? " checkInside" : " checkOutside";
  std::cout << name << " test: passed" << std::endl;
  return true;
}

static bool testLocalMethod ()
{
  bool pass = true;
  static const int cdim = 2;
  static const int mydim = 2;
  typedef double ctype;

  typedef Dune::MultiLinearGeometry< ctype, mydim, cdim, Dune::MultiLinearGeometryTraits< ctype > > Geometry;

  {
    Dune::GeometryType type = Dune::GeometryTypes::simplex( cdim );
    auto refElement = Dune::referenceElement< ctype, cdim >( type );
    const int numCorners = refElement.size( mydim );
    std::vector< Dune::FieldVector< ctype, cdim > > corners( numCorners );
    corners[ 0 ] = {{ 0.5, 0 }};
    corners[ 1 ] = {{ 0.5, 0.5 }};
    corners[ 2 ] = {{ 0  , 0.5 }};

    Geometry geometry( refElement, corners );
    // point inside of geometry
    Dune::FieldVector< ctype, cdim > inside( 0.4 );
    pass &= testCheckInside( refElement, geometry, inside, true );

    // point outside of geometry
    Dune::FieldVector< ctype, cdim > outside( -0.1 );
    pass &= testCheckInside( refElement, geometry, outside, false );

    // point outside of geometry
    Dune::FieldVector< ctype, cdim > outside2( 0 );
    pass &= testCheckInside( refElement, geometry, outside2, false );
  }

  {
    Dune::GeometryType type = Dune::GeometryTypes::cube( cdim );
    auto refElement = Dune::referenceElement< ctype, cdim >( type );
    const int numCorners = refElement.size( mydim );
    std::vector< Dune::FieldVector< ctype, cdim > > corners( numCorners );
    corners[ 0 ] = {{ 0.0, 0.0 }};
    corners[ 1 ] = {{ 0.5, 0.5 }};
    corners[ 2 ] = {{ 0.0, 1.0 }};
    corners[ 3 ] = {{ 0.5, 1.0 }};

    Geometry geometry( refElement, corners );
    // point inside
    Dune::FieldVector< ctype, cdim > inside( {{ 0.25, 0.75 }});
    pass &= testCheckInside( refElement, geometry, inside, true );
    // point outside of geometry
    Dune::FieldVector< ctype, cdim > outside( {{ 1, 0}} );
    pass &= testCheckInside( refElement, geometry, outside, false );
    // another point outside
    Dune::FieldVector< ctype, cdim > outside2( 4 );
    pass &= testCheckInside( refElement, geometry, outside2, false );
  }
  return pass;
}


template< class ctype, int mydim, int cdim, class Traits >
static bool testMultiLinearGeometry ( Dune::Transitional::ReferenceElement< ctype, Dune::Dim<mydim> > refElement,
                                      const Dune::FieldMatrix< ctype, mydim, mydim > &A,
                                      const Dune::FieldMatrix< ctype, cdim, cdim > &B,
                                      const Traits & /* traits */ )
{
  bool pass = true;

  typedef Dune::MultiLinearGeometry< ctype, mydim, cdim, Traits > Geometry;

  const Dune::FieldVector< ctype, mydim > &localCenter = refElement.position( 0, 0 );
  const ctype epsilon = ctype( 1e5 )*std::numeric_limits< ctype >::epsilon();

  const ctype detA = A.determinant();
  assert( std::abs( std::abs( B.determinant() ) - ctype( 1 ) ) <= epsilon );

  const int numCorners = refElement.size( mydim );
  std::vector< Dune::FieldVector< ctype, cdim > > corners( numCorners );
  for( int i = 0; i < numCorners; ++i )
    corners[ i ] = map( A, B, refElement.position( i, mydim ) );

  Geometry geometry( refElement, corners );

  if( !geometry.affine() )
  {
    std::cerr << "Error: Affine returns false for an affine geometry." << std::endl;
    pass = false;
  }

  const ctype integrationElement = geometry.integrationElement( localCenter );
  if( std::abs( integrationElement - std::abs( detA ) ) > epsilon )
  {
    std::cerr << "Error: Wrong integration element (" << integrationElement
              << ", should be " << std::abs( detA )
              << ")." << std::endl;
    pass = false;
  }
  const ctype volume = geometry.volume();
  if( std::abs( volume - refElement.volume()*std::abs( detA ) ) > epsilon )
  {
    std::cerr << "Error: Wrong volume (" << volume
              << ", should be " << (refElement.volume()*std::abs( detA ))
              << ")." << std::endl;
    pass = false;
  }

  const Dune::FieldVector< ctype, cdim > center = geometry.center();
  if( (center - map( A, B, refElement.position( 0, 0 ) )).two_norm() > epsilon )
  {
    std::cerr << "Error: wrong barycenter (" << center << ")." << std::endl;
    pass = false;
  }

  for (int c = 0; c < numCorners; ++c) {
    Dune::FieldVector<ctype, mydim> local(refElement.position(c, mydim));
    Dune::FieldVector<ctype, cdim> global(geometry.global(local));
    Dune::FieldVector<ctype, mydim> local2(geometry.local(global));
    if (local2 != local2) {
      std::cerr << "Error: at corner " << c << " local returned NaN: "
                << local2 << std::endl;
      pass = false;
    }
    if ((local - local2).two_norm() > epsilon) {
      std::cerr << "Error: at corner " << c << " local returned wrong value: "
                << local2 << " (expected: " << local << ")" << std::endl;
      pass = false;
    }
  }

  const Dune::FieldMatrix< ctype, mydim, cdim > JT = geometry.jacobianTransposed( localCenter );
  for( int i = 0; i < mydim; ++i )
  {
    Dune::FieldVector< ctype, mydim > e( ctype( 0 ) );
    e[ i ] = ctype( 1 );
    const Dune::FieldVector< ctype, cdim > t = map( A, B, e );
    if (JT[i] != JT[i])
    {
      std::cerr << "Error: jacobianTransposed[" << i << "] (" << JT[i]
                << ") has NaN entry." << std::endl;
      pass = false;
    }
    if( (t - JT[ i ]).two_norm() > epsilon )
    {
      std::cerr << "Error: wrong jacobianTransposed[ " << i << " ] (" << JT[ i ]
                << ", should be " << t << ")." << std::endl;
      pass = false;
    }
  }

  for (int c = 0; c < numCorners; ++c) {
    const Dune::FieldMatrix< ctype, mydim, cdim > cornerJT = geometry.jacobianTransposed(refElement.position(c, mydim));
    for( int i = 0; i < mydim; ++i ) {
      Dune::FieldVector< ctype, mydim > e( ctype( 0 ) );
      e[ i ] = ctype( 1 );
      const Dune::FieldVector< ctype, cdim > t = map( A, B, e );
      if (cornerJT[i] != cornerJT[i])
      {
        std::cerr << "Error: at corner " << c << ": jacobianTransposed["
                  << i << "] (" << cornerJT[i] << ") has NaN entry." << std::endl;
        pass = false;
      }
      if( (t - cornerJT[ i ]).two_norm() > epsilon )
      {
        std::cerr << "Error: at corner " << c << ": wrong jacobianTransposed[ "
                  << i << " ] (" << cornerJT[ i ] << ", should be " << t << ")."
                  << std::endl;
        pass = false;
      }
    }
  }

  pass &= checkGeometry( geometry );

  return pass;
}


template< class ctype, int mydim, int cdim, class Traits >
static bool testMultiLinearGeometry ( Dune::GeometryType gt,
                                      const Traits &traits )
{
  auto refElement = Dune::referenceElement< ctype, mydim >( gt );

  Dune::FieldMatrix< ctype, mydim, mydim > A;
  Dune::FieldMatrix< ctype, cdim, cdim > B;

  std::cout << "Checking geometry (topologyId = " << gt.id() << ", mydim = " << mydim << ", cdim = " << cdim << ")" << std::endl;
  std::cout << " reference mapping: ";
  A = ctype( 0 );
  for( int i = 0; i < mydim; ++i )
    A[ i ][ i ] = ctype( 1 );
  B = ctype( 0 );
  for( int i = 0; i < cdim; ++i )
    B[ i ][ i ] = ctype( 1 );
  const bool passId = testMultiLinearGeometry( refElement, A, B, traits );
  std::cout << (passId ? "passed" : "failed");

  std::cout << ", scaled reference mapping: ";
  A = ctype( 0 );
  for( int i = 0; i < mydim; ++i )
    A[ i ][ i ] = ctype( 42 );
  B = ctype( 0 );
  for( int i = 0; i < cdim; ++i )
    B[ i ][ i ] = ctype( 1 );
  const bool passScaledId = testMultiLinearGeometry( refElement, A, B, traits );
  std::cout << (passScaledId ? "passed" : "failed") << std::endl;

  return passId && passScaledId;
}

template<class ctype, class Traits>
static bool testNonLinearGeometry(const Traits & /* traits */)
{
  const unsigned dim = 2;
  typedef Dune::FieldVector<ctype,dim> Vector;
  typedef Dune::MultiLinearGeometry<ctype,dim,dim,Traits> Geometry;
  const ctype epsilon(ctype(16) * std::numeric_limits<ctype>::epsilon());

  bool pass(true);
  std::cout << "Checking geometry (non-linear, quadrilateral): ";

  auto reference = Dune::referenceElement<ctype,dim>(Dune::GeometryTypes::quadrilateral);

  std::vector<Vector> corners = {{0,0},
                                 {2,0},
                                 {0,1},
                                 {1,1}};

  const Geometry geometry(reference, corners);

  /* Test global() */
  for (std::size_t c = 0; c < corners.size(); ++c) {
    const Vector& local(reference.position(c, dim));
    const Vector global(geometry.global(local));
    if (global != global) {
      std::cerr << "global failed at corner " << c << ": returned NaN: "
                << global << std::endl;
      pass = false;
    }
    if ((global - corners[c]).two_norm() > epsilon) {
      std::cerr << "global failed at corner " << c << ": got " << global
                << ", but expected " << corners[c] << std::endl;
      pass = false;
    }
  }

  /* Test global() outside reference element */
  {
    Vector local  = {-1, 0};
    Vector global = {-2, 0};
    const Vector global2(geometry.global(local));
    if (global2 != global2) {
      std::cerr << "global failed outside reference element: returned NaN: "
                << global2 << std::endl;
      pass = false;
    }
    if ((global - global2).two_norm() > epsilon) {
      std::cerr << "global failed outside reference element: got " << global2
                << ", but expected " << global << std::endl;
      pass = false;
    }
  }

  /* Test local() */
  for (std::size_t c = 0; c < corners.size(); ++c) {
    const Vector& local(reference.position(c, dim));
    const Vector local2(geometry.local(corners[c]));
    if (local2 != local2) {
      std::cerr << "local failed at corner " << c << ": returned NaN: "
                << local2 << std::endl;
      pass = false;
    }
    if ((local - local2).two_norm() > epsilon) {
      std::cerr << "local failed at corner " << c << ": got " << local2
                << ", but expected " << local << std::endl;
      pass = false;
    }
  }

  /* Test local() outside reference element */
  {
    Vector global = {-2, 0};
    Vector local  = {-1, 0};
    const Vector local2(geometry.local(global));
    if (local2 != local2) {
      std::cerr << "local failed outside reference element: returned NaN: "
                << local2 << std::endl;
      pass = false;
    }
    if ((local - local2).two_norm() > epsilon) {
      std::cerr << "local failed outside reference element: got " << local2
                << ", but expected " << local << std::endl;
      pass = false;
    }
  }

  std::cout << (pass ? "passed" : "failed") << std::endl;
  return pass;
}

template< class ctype, class Traits >
static bool testMultiLinearGeometry ( const Traits& traits )
{
  bool pass = true;

  Dune::GeometryType cube0d = Dune::GeometryTypes::cube(0);
  Dune::GeometryType cube1d = Dune::GeometryTypes::cube(1);
  Dune::GeometryType cube2d = Dune::GeometryTypes::cube(2);
  Dune::GeometryType cube3d = Dune::GeometryTypes::cube(3);
  Dune::GeometryType cube4d = Dune::GeometryTypes::cube(4);

  Dune::GeometryType simplex0d = Dune::GeometryTypes::simplex(0);
  Dune::GeometryType simplex1d = Dune::GeometryTypes::simplex(1);
  Dune::GeometryType simplex2d = Dune::GeometryTypes::simplex(2);
  Dune::GeometryType simplex3d = Dune::GeometryTypes::simplex(3);
  Dune::GeometryType simplex4d = Dune::GeometryTypes::simplex(4);

  Dune::GeometryType prism3d = Dune::GeometryTypes::prism;
  Dune::GeometryType pyramid3d = Dune::GeometryTypes::pyramid;

  pass &= testMultiLinearGeometry< ctype, 0, 0 >( simplex0d, traits );
  pass &= testMultiLinearGeometry< ctype, 0, 1 >( simplex0d, traits );
  pass &= testMultiLinearGeometry< ctype, 0, 2 >( simplex0d, traits );
  pass &= testMultiLinearGeometry< ctype, 0, 3 >( simplex0d, traits );
  pass &= testMultiLinearGeometry< ctype, 0, 4 >( simplex0d, traits );

  pass &= testMultiLinearGeometry< ctype, 0, 0 >( cube0d, traits );
  pass &= testMultiLinearGeometry< ctype, 0, 1 >( cube0d, traits );
  pass &= testMultiLinearGeometry< ctype, 0, 2 >( cube0d, traits );
  pass &= testMultiLinearGeometry< ctype, 0, 3 >( cube0d, traits );
  pass &= testMultiLinearGeometry< ctype, 0, 4 >( cube0d, traits );

  pass &= testMultiLinearGeometry< ctype, 1, 1 >( simplex1d, traits );
  pass &= testMultiLinearGeometry< ctype, 1, 2 >( simplex1d, traits );
  pass &= testMultiLinearGeometry< ctype, 1, 3 >( simplex1d, traits );
  pass &= testMultiLinearGeometry< ctype, 1, 4 >( simplex1d, traits );

  pass &= testMultiLinearGeometry< ctype, 1, 3 >( cube1d, traits );
  pass &= testMultiLinearGeometry< ctype, 1, 1 >( cube1d, traits );
  pass &= testMultiLinearGeometry< ctype, 1, 2 >( cube1d, traits );
  pass &= testMultiLinearGeometry< ctype, 1, 4 >( cube1d, traits );

  pass &= testMultiLinearGeometry< ctype, 2, 2 >( simplex2d, traits );
  pass &= testMultiLinearGeometry< ctype, 2, 3 >( simplex2d, traits );
  pass &= testMultiLinearGeometry< ctype, 2, 4 >( simplex2d, traits );

  pass &= testMultiLinearGeometry< ctype, 2, 2 >( cube2d, traits );
  pass &= testMultiLinearGeometry< ctype, 2, 3 >( cube2d, traits );
  pass &= testMultiLinearGeometry< ctype, 2, 4 >( cube2d, traits );

  pass &= testMultiLinearGeometry< ctype, 3, 3 >( simplex3d, traits );
  pass &= testMultiLinearGeometry< ctype, 3, 4 >( simplex3d, traits );

  pass &= testMultiLinearGeometry< ctype, 3, 3 >( pyramid3d, traits );
  pass &= testMultiLinearGeometry< ctype, 3, 4 >( pyramid3d, traits );

  pass &= testMultiLinearGeometry< ctype, 3, 3 >( prism3d, traits );
  pass &= testMultiLinearGeometry< ctype, 3, 4 >( prism3d, traits );

  pass &= testMultiLinearGeometry< ctype, 3, 3 >( cube3d, traits );
  pass &= testMultiLinearGeometry< ctype, 3, 4 >( cube3d, traits );

  pass &= testMultiLinearGeometry< ctype, 4, 4 >( simplex4d, traits );
  pass &= testMultiLinearGeometry< ctype, 4, 5 >( simplex4d, traits );

  pass &= testMultiLinearGeometry< ctype, 4, 4 >( cube4d, traits );
  pass &= testMultiLinearGeometry< ctype, 4, 5 >( cube4d, traits );

  pass &= testNonLinearGeometry<ctype>( traits );

  return pass;
}

int main ( int /* argc */, char ** /* argv */)
{
  bool pass = true;

  std::cout << ">>> Checking ctype = double" << std::endl;
  pass &= testMultiLinearGeometry< double >
    ( Dune::MultiLinearGeometryTraits< double >{} );

  std::cout << ">>> Checking ctype = double with reference_wrapped corner "
            << "storage" << std::endl;
  pass &= testMultiLinearGeometry< double >
    ( ReferenceWrapperGeometryTraits< double >{} );

  // std::cout << ">>> Checking ctype = float" << std::endl;
  // pass &= testMultiLinearGeometry< float >
  //   ( Dune::MultiLinearGeometryTraits< float >{} );

  // std::cout << ">>> Checking ctype = float with reference_wrapped corner "
  //           << "storage" << std::endl;
  // pass &= testMultiLinearGeometry< float >
  //   ( ReferenceWrapperGeometryTraits< float >{} );

  pass &= testLocalMethod();

  return (pass ? 0 : 1);
}
