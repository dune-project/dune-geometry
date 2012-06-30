// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_CHECK_GEOMETRY_HH
#define DUNE_CHECK_GEOMETRY_HH

#include <limits>

#include <dune/common/typetraits.hh>

#include <dune/geometry/quadraturerules/gaussquadrature.hh>

namespace Dune
{

  /** \brief Static and dynamic checks for all features of a Geometry
   *
   * This excludes anything related to being part of a grid.
   *
   * \tparam TestGeometry The type of the geometry to be tested
   *
   * \param geometry The TestGeometry object to be tested
   *
   * \returns true if check passed
   */
  template <class TestGeometry>
  bool checkGeometry ( const TestGeometry& geometry )
  {
    bool pass = true;

    ////////////////////////////////////////////////////////////////
    // Extract all static information
    ////////////////////////////////////////////////////////////////

    // dimension of the corresponding reference element
    static const int mydim = TestGeometry::mydimension;

    // dimension of the world space
    static const int coorddim = TestGeometry::coorddimension;

    // type used for coordinate coefficients
    typedef typename TestGeometry::ctype ctype;

    // vector type used for points in the domain
    typedef typename TestGeometry::LocalCoordinate LocalCoordinate;

    // vector type used for image points
    typedef typename TestGeometry::GlobalCoordinate GlobalCoordinate;

    ////////////////////////////////////////////////////////////////////////

    GeometryType type = geometry.type();

    const GenericReferenceElement< ctype, mydim > &refElement = GenericReferenceElements< ctype, mydim >::general(type);

    // Test whether the return value of the method 'center' corresponds to the center of the
    // reference element.  That is the current definition of the method.
    const FieldVector<ctype, coorddim> center = geometry.global(refElement.position(0,0));
    if( std::abs( (geometry.center() - center).two_norm() ) > 1e-8 )
      DUNE_THROW(Exception, "center() is not consistent with global(refElem.position(0,0)).");

    ////////////////////////////////////////////////////////////////////////
    // Test whether the number and placement of the corners is consistent
    // with the corners of the corresponding reference element.
    ////////////////////////////////////////////////////////////////////////

    if( refElement.size( mydim ) == geometry.corners() )
    {
      for( int i = 0; i < geometry.corners(); ++i )
      {
        if( (geometry.corner( i ) - geometry.global( refElement.position( i, mydim ) )).two_norm() > 1e-8 ) {
          std::cerr << "Error: Methods corner and global are inconsistent." << std::endl;
          pass = false;
        }
      }
    }
    else {
      std::cerr << "Error: Incorrect number of corners (" << geometry.corners() << ", should be " << refElement.size( mydim ) << ")." << std::endl;
      pass = false;
    }

    ///////////////////////////////////////////////////////////////////////////////
    // Use a quadrature rule as a set of test points and loop over them
    ///////////////////////////////////////////////////////////////////////////////
    typedef Dune::GenericGeometry::GaussPoints< ctype > Points;
    typedef Dune::GenericGeometry::GenericQuadratureFactory< mydim, ctype, Points > QuadratureFactory;
    const typename QuadratureFactory::Object &quadrature = *QuadratureFactory::create( geometry.type(), 2 );
    for( size_t i = 0; i < quadrature.size(); ++i )
    {
      const typename TestGeometry::LocalCoordinate &x = quadrature[ i ].position();

      // Test whether the methods 'local' and 'global' are inverse to each other
      if ( (x - geometry.local( geometry.global( x ) )).two_norm() > 1e-8 ) {
        std::cerr << "Error: global and local are not inverse to each other." << std::endl;
        pass = false;
      }

      // Test whether the methods 'jacobianTransposed' and 'jacobianInverseTransposed'
      // return matrices that are inverse to each other.
      const FieldMatrix< ctype, mydim, coorddim > &jt = geometry.jacobianTransposed( x );
      const FieldMatrix< ctype, coorddim, mydim > &jit = geometry.jacobianInverseTransposed( x );

      FieldMatrix< ctype, mydim, mydim > id;
      FMatrixHelp::multMatrix( jt, jit, id );
      bool isId = true;
      for( int j = 0; j < mydim; ++j )
        for( int k = 0; k < mydim; ++k )
          isId &= (std::abs( id[ j ][ k ] - (j == k ? 1 : 0) ) < 1e-8);
      if( !isId)
      {
        std::cerr << "Error: jacobianTransposed and jacobianInverseTransposed are not inverse to each other." << std::endl;
        std::cout << "       id != [ ";
        for( int j = 0; j < mydim; ++j )
          std::cout << (j > 0 ? " | " : "") << id[ j ];
        std::cout << " ]" << std::endl;
        pass = false;
      }

      // Test whether integrationElement returns something nonnegative
      if( geometry.integrationElement( x ) < 0 ) {
        std::cerr << "Error: Negative integrationElement found." << std::endl;
        pass = false;
      }

      FieldMatrix< ctype, mydim, mydim > jtj( 0 );
      for( int i = 0; i < mydim; ++i )
        for( int j = 0; j < mydim; ++j )
          for( int k = 0; k < coorddim; ++k )
            jtj[ i ][ j ] += jt[ i ][ k ] * jt[ j ][ k ];
      if( std::abs( std::sqrt( jtj.determinant() ) - geometry.integrationElement( x ) ) > 1e-8 ) {
        std::cerr << "Error: integrationElement is not consistent with jacobianTransposed." << std::endl;
        pass = false;
      }
      if (geometry.affine())
        if( std::abs( geometry.volume() - refElement.volume()*geometry.integrationElement( x ) ) > 1e-8 ) {
          std::cerr << "Error: volume is not consistent with jacobianTransposed." << std::endl;
          pass = false;
        }
    }
    QuadratureFactory::release( &quadrature );

    return pass;
  }

}

#endif // #ifndef DUNE_CHECK_GEOMETRY_HH
