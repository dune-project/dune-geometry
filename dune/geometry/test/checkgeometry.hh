// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_CHECK_GEOMETRY_HH
#define DUNE_CHECK_GEOMETRY_HH

#include <limits>

#include <dune/common/typetraits.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

namespace Dune
{

  /**
   * \brief Static and dynamic checks for all features of a Geometry
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
    using std::sqrt;
    using std::abs;
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
    [[maybe_unused]] typedef typename TestGeometry::LocalCoordinate LocalCoordinate;

    // vector type used for image points
    [[maybe_unused]] typedef typename TestGeometry::GlobalCoordinate GlobalCoordinate;

    // Matrix-like type for the return value of the jacobianTransposed method
    typedef typename TestGeometry::JacobianTransposed JacobianTransposed;

    // Matrix-like type for the return value of the jacobianInverseTransposed method
    typedef typename TestGeometry::JacobianInverseTransposed JacobianInverseTransposed;

    const ctype tolerance = std::sqrt( std::numeric_limits< ctype >::epsilon() );

    ////////////////////////////////////////////////////////////////////////

    const int corners = geometry.corners();
    if( corners == 0 )
      DUNE_THROW( Exception, "A geometry must have at least one corner." );

    GlobalCoordinate cornerAvg ( 0 );
    for( int i = 0; i < corners; ++i )
      cornerAvg += geometry.corner( i );
    cornerAvg /= ctype( corners );

    const GlobalCoordinate center = geometry.center();

    if( corners > 1 )
    {
      // if we have more than one corner, no corner may coincide with their average or the center
      for( int i = 0; i < corners; ++i )
      {
        const GlobalCoordinate corner = geometry.corner( i );
        if( (corner - center).two_norm() <= tolerance )
        {
          std::cerr << "Error: Geometry has multiple corners, but one corner coincides with the center." << std::endl;
          pass = false;
        }

        if( (corner - cornerAvg).two_norm() <= tolerance )
        {
          std::cerr << "Error: Geometry has multiple corners, but one corner coincides with their average." << std::endl;
          pass = false;
        }
      }
    }
    else
    {
      // the single corner must coincide with the center
      if( (center - cornerAvg).two_norm() > tolerance )
      {
        std::cerr << "Error: Geometry has a single corner (" << cornerAvg << "), but it does not coincide with the center (" << center << ")." << std::endl;
        pass = false;
      }
    }

    const ctype volume = geometry.volume();
    if( volume < tolerance )
    {
      std::cerr << "Error: Geometry has nearly vanishing volume (" << volume << ")" << std::endl;
      pass = false;
    }

    ////////////////////////////////////////////////////////////////////////

    const GeometryType type = geometry.type();
    if( type.isNone() )
      return pass;

    // make sure the reference element type lookup works
    ReferenceElement< TestGeometry > refElement = referenceElement( geometry );

    // Test whether the return value of the method 'center' corresponds to the center of the
    // reference element.  That is the current definition of the method.
    if( (center - geometry.global( refElement.position( 0, 0 ) )).two_norm() > tolerance )
      DUNE_THROW( Exception, "center() is not consistent with global(refElem.position(0,0))." );

    ////////////////////////////////////////////////////////////////////////
    // Test whether the number and placement of the corners is consistent
    // with the corners of the corresponding reference element.
    ////////////////////////////////////////////////////////////////////////

    if( refElement.size( mydim ) == corners )
    {
      for( int i = 0; i < geometry.corners(); ++i )
      {
        if( (geometry.corner( i ) - geometry.global( refElement.position( i, mydim ) )).two_norm() > tolerance )
        {
          std::cerr << "Error: Methods corner and global are inconsistent." << std::endl;
          pass = false;
        }
      }
    }
    else
    {
      std::cerr << "Error: Incorrect number of corners (" << geometry.corners() << ", should be " << refElement.size( mydim ) << ")." << std::endl;
      pass = false;
    }

    ///////////////////////////////////////////////////////////////////////////////
    // Use a quadrature rule as a set of test points and loop over them
    ///////////////////////////////////////////////////////////////////////////////
    const Dune::QuadratureRule<ctype, mydim> & quadrature = Dune::QuadratureRules<ctype, mydim>::rule(geometry.type(), 2);
    for (const auto& ip : quadrature)
    {
      const typename TestGeometry::LocalCoordinate &x = ip.position();

      // Test whether the methods 'local' and 'global' are inverse to each other
      if ( (x - geometry.local( geometry.global( x ) )).two_norm() > tolerance ) {
        std::cerr << "Error: global and local are not inverse to each other." << std::endl;
        pass = false;
      }

      // Test whether the methods 'jacobianTransposed' and 'jacobianInverseTransposed'
      // return matrices that are inverse to each other.
      const JacobianTransposed &jt = geometry.jacobianTransposed( x );
      const JacobianInverseTransposed &jit = geometry.jacobianInverseTransposed( x );

      // Transform to FieldMatrix, so we can have coefficent access and other goodies
      // We need some black magic for the transformation, because there is no
      // official easy way yet.
      // The following code does the transformation by multiplying jt and jit from
      // the right by identity matrices.  That way, only the mv method is used.
      FieldMatrix< ctype, mydim, coorddim > jtAsFieldMatrix;
      for (int j=0; j<coorddim; j++) {

        FieldVector<ctype,coorddim> idColumn(0);
        idColumn[j] = 1;

        FieldVector<ctype,mydim> column;
        jt.mv(idColumn,column);

        for (int k=0; k<mydim; k++)
          jtAsFieldMatrix[k][j] = column[k];

      }

      FieldMatrix< ctype, coorddim, mydim > jitAsFieldMatrix;
      for (int j=0; j<mydim; j++) {

        FieldVector<ctype,mydim> idColumn(0);
        idColumn[j] = 1;

        FieldVector<ctype,coorddim> column;
        jit.mv(idColumn,column);

        for (int k=0; k<coorddim; k++)
          jitAsFieldMatrix[k][j] = column[k];

      }


      FieldMatrix< ctype, mydim, mydim > id;
      FMatrixHelp::multMatrix( jtAsFieldMatrix, jitAsFieldMatrix, id );
      bool isId = true;
      for( int j = 0; j < mydim; ++j )
        for( int k = 0; k < mydim; ++k )
          isId &= (std::abs( id[ j ][ k ] - (j == k ? 1 : 0) ) < tolerance);
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
            jtj[ i ][ j ] += jtAsFieldMatrix[ i ][ k ] * jtAsFieldMatrix[ j ][ k ];

      if( abs( sqrt( jtj.determinant() ) - geometry.integrationElement( x ) ) > tolerance ) {
        std::cerr << "Error: integrationElement is not consistent with jacobianTransposed." << std::endl;
        pass = false;
      }
      if (geometry.affine())
        if( abs( geometry.volume() - refElement.volume()*geometry.integrationElement( x ) ) > tolerance ) {
          std::cerr << "Error: volume is not consistent with jacobianTransposed." << std::endl;
          pass = false;
        }
    }

    return pass;
  }

}

#endif // #ifndef DUNE_CHECK_GEOMETRY_HH
