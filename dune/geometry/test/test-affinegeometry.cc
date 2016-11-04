// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <vector>

#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/geometry/test/checkgeometry.hh>


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


template< class ctype, int mydim, int cdim >
static bool testAffineGeometry ( const Dune::ReferenceElement< ctype, mydim > &refElement,
                                 const Dune::FieldMatrix< ctype, mydim, mydim > &A,
                                 const Dune::FieldMatrix< ctype, cdim, cdim > &B )
{
  bool pass = true;

  typedef Dune::AffineGeometry< ctype, mydim, cdim > Geometry;

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

  const Dune::FieldMatrix< ctype, mydim, cdim > JT = geometry.jacobianTransposed( localCenter );
  for( int i = 0; i < mydim; ++i )
  {
    Dune::FieldVector< ctype, mydim > e( ctype( 0 ) );
    e[ i ] = ctype( 1 );
    const Dune::FieldVector< ctype, cdim > t = map( A, B, e );
    if( (t - JT[ i ]).two_norm() > epsilon )
    {
      std::cerr << "Error: wrong jacobianTransposed[ " << i << " ] (" << JT[ i ]
                << ", should be " << t << ")." << std::endl;
      pass = false;
    }

  }

  pass &= checkGeometry( geometry );

  return pass;
}


template< class ctype, int mydim, int cdim >
static bool testAffineGeometry ( Dune::GeometryType gt )
{
  const Dune::ReferenceElement< ctype, mydim > &refElement = Dune::ReferenceElements< ctype, mydim >::general( gt );

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
  const bool passId = testAffineGeometry( refElement, A, B );
  std::cout << (passId ? "passed" : "failed");

  std::cout << ", scaled reference mapping: ";
  A = ctype( 0 );
  for( int i = 0; i < mydim; ++i )
    A[ i ][ i ] = ctype( 42 );
  B = ctype( 0 );
  for( int i = 0; i < cdim; ++i )
    B[ i ][ i ] = ctype( 1 );
  const bool passScaledId = testAffineGeometry( refElement, A, B );
  std::cout << (passScaledId ? "passed" : "failed") << std::endl;

  return passId && passScaledId;
}


template< class ctype >
static bool testAffineGeometry ()
{
  bool pass = true;

  Dune::GeometryType cube0d( Dune::Impl::CubeTopology< 0 >::type::id, 0 );
  Dune::GeometryType cube1d( Dune::Impl::CubeTopology< 1 >::type::id, 1 );
  Dune::GeometryType cube2d( Dune::Impl::CubeTopology< 2 >::type::id, 2 );
  Dune::GeometryType cube3d( Dune::Impl::CubeTopology< 3 >::type::id, 3 );
  Dune::GeometryType cube4d( Dune::Impl::CubeTopology< 4 >::type::id, 4 );

  Dune::GeometryType simplex0d( Dune::Impl::SimplexTopology< 0 >::type::id, 0 );
  Dune::GeometryType simplex1d( Dune::Impl::SimplexTopology< 1 >::type::id, 1 );
  Dune::GeometryType simplex2d( Dune::Impl::SimplexTopology< 2 >::type::id, 2 );
  Dune::GeometryType simplex3d( Dune::Impl::SimplexTopology< 3 >::type::id, 3 );
  Dune::GeometryType simplex4d( Dune::Impl::SimplexTopology< 4 >::type::id, 4 );

  Dune::GeometryType prism3d( Dune::Impl::PrismTopology< 3 >::type::id, 3 );
  Dune::GeometryType pyramid3d( Dune::Impl::PyramidTopology< 3 >::type::id, 3 );

  pass &= testAffineGeometry< ctype, 0, 0 >( simplex0d );
  pass &= testAffineGeometry< ctype, 0, 1 >( simplex0d );
  pass &= testAffineGeometry< ctype, 0, 2 >( simplex0d );
  pass &= testAffineGeometry< ctype, 0, 3 >( simplex0d );
  pass &= testAffineGeometry< ctype, 0, 4 >( simplex0d );

  pass &= testAffineGeometry< ctype, 0, 0 >( cube0d );
  pass &= testAffineGeometry< ctype, 0, 1 >( cube0d );
  pass &= testAffineGeometry< ctype, 0, 2 >( cube0d );
  pass &= testAffineGeometry< ctype, 0, 3 >( cube0d );
  pass &= testAffineGeometry< ctype, 0, 4 >( cube0d );

  pass &= testAffineGeometry< ctype, 1, 1 >( simplex1d );
  pass &= testAffineGeometry< ctype, 1, 2 >( simplex1d );
  pass &= testAffineGeometry< ctype, 1, 3 >( simplex1d );
  pass &= testAffineGeometry< ctype, 1, 4 >( simplex1d );

  pass &= testAffineGeometry< ctype, 1, 3 >( cube1d );
  pass &= testAffineGeometry< ctype, 1, 1 >( cube1d );
  pass &= testAffineGeometry< ctype, 1, 2 >( cube1d );
  pass &= testAffineGeometry< ctype, 1, 4 >( cube1d );

  pass &= testAffineGeometry< ctype, 2, 2 >( simplex2d );
  pass &= testAffineGeometry< ctype, 2, 3 >( simplex2d );
  pass &= testAffineGeometry< ctype, 2, 4 >( simplex2d );

  pass &= testAffineGeometry< ctype, 2, 2 >( cube2d );
  pass &= testAffineGeometry< ctype, 2, 3 >( cube2d );
  pass &= testAffineGeometry< ctype, 2, 4 >( cube2d );

  pass &= testAffineGeometry< ctype, 3, 3 >( simplex3d );
  pass &= testAffineGeometry< ctype, 3, 4 >( simplex3d );

  /** \bug These tests currently fail. */
  //   pass &= testAffineGeometry< ctype, 3, 3 >( pyramid3d );
  //   pass &= testAffineGeometry< ctype, 3, 4 >( pyramid3d );

  pass &= testAffineGeometry< ctype, 3, 3 >( prism3d );
  pass &= testAffineGeometry< ctype, 3, 4 >( prism3d );

  /** \bug These tests currently fail. */
  //   pass &= testAffineGeometry< ctype, 3, 3 >( cube3d );
  //   pass &= testAffineGeometry< ctype, 3, 4 >( cube3d );

  pass &= testAffineGeometry< ctype, 4, 4 >( simplex4d );
  pass &= testAffineGeometry< ctype, 4, 5 >( simplex4d );

  /** \bug These tests currently fail. */
  //   pass &= testAffineGeometry< ctype, 4, 4 >( cube4d );
  //   pass &= testAffineGeometry< ctype, 4, 5 >( cube4d );

  return pass;
}

int main ( int argc, char **argv )
{
  bool pass = true;

  std::cout << ">>> Checking ctype = double" << std::endl;
  pass &= testAffineGeometry< double >();
  //std::cout << ">>> Checking ctype = float" << std::endl;
  //pass &= testAffineGeometry< float >();

  return (pass ? 0 : 1);
}
