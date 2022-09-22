// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

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

  namespace Impl
  {

    // Helper function to convert a matrix to a FieldMatrix
    // using only the mv method.
    template<class ctype, int rows, int cols, class M>
    Dune::FieldMatrix<ctype, rows, cols> toFieldMatrix(const M& m){
      Dune::FieldMatrix<ctype, rows, cols> result;
      for (int j=0; j<cols; j++) {
        FieldVector<ctype,cols> idColumn(0);
        idColumn[j] = 1;
        FieldVector<ctype,rows> column;
        m.mv(idColumn,column);
        for (int k=0; k<rows; k++)
          result[k][j] = column[k];
      }
      return result;
    }

    // Helper function to check if FieldMatrix is an identity matrix
    template<class ctype, int rows, int cols>
    bool isIdentity(const Dune::FieldMatrix<ctype, rows, cols>& m, double tolerance){
      if (rows != cols)
        return false;
      for(int i=0; i<rows; ++i)
        for(int j=0; j<rows; ++j)
          if (std::abs(m[i][j] - (i == j ? 1 : 0)) >= tolerance)
            return false;
      return true;
    }

  }

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

    // Matrix-like type for the return value of the jacobian method
    typedef typename TestGeometry::Jacobian Jacobian;

    // Matrix-like type for the return value of the jacobianInverse method
    typedef typename TestGeometry::JacobianInverse JacobianInverse;

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

      const JacobianTransposed &Jt = geometry.jacobianTransposed( x );
      const JacobianInverseTransposed &Jit = geometry.jacobianInverseTransposed( x );
      const Jacobian &J = geometry.jacobian( x );
      const JacobianInverse &Ji = geometry.jacobianInverse( x );

      // Transform to FieldMatrix, so we can have coefficient access and other goodies
      auto JtAsFieldMatrix = Impl::toFieldMatrix< ctype, mydim, coorddim >(Jt);
      auto JitAsFieldMatrix = Impl::toFieldMatrix< ctype, coorddim, mydim >(Jit);
      auto JAsFieldMatrix = Impl::toFieldMatrix< ctype, coorddim, mydim >(J);
      auto JiAsFieldMatrix = Impl::toFieldMatrix< ctype, mydim, coorddim >(Ji);


      // Test whether the methods 'jacobianTransposed' and 'jacobianInverseTransposed'
      // return matrices that are inverse to each other.
      {
        FieldMatrix< ctype, mydim, mydim > id = JtAsFieldMatrix * JitAsFieldMatrix;
        if(not Impl::isIdentity(id, tolerance))
        {
          std::cerr << "Error: jacobianTransposed and jacobianInverseTransposed are not inverse to each other." << std::endl;
          std::cout << "       id != [ ";
          for( int j = 0; j < mydim; ++j )
            std::cout << (j > 0 ? " | " : "") << id[ j ];
          std::cout << " ]" << std::endl;
          pass = false;
        }
      }

      // Test whether the methods 'jacobian' and 'jacobianInverse'
      // return matrices that are inverse to each other.
      {
        FieldMatrix< ctype, mydim, mydim > id = JiAsFieldMatrix * JAsFieldMatrix;
        if(not Impl::isIdentity(id, tolerance))
        {
          std::cerr << "Error: jacobian and jacobianInverse are not inverse to each other." << std::endl;
          std::cout << "       id != [ ";
          for( int j = 0; j < mydim; ++j )
            std::cout << (j > 0 ? " | " : "") << id[ j ];
          std::cout << " ]" << std::endl;
          pass = false;
        }
      }

      // Test whether the methods 'jacobianTransposed' and 'jacobianInverseTransposed'
      // are the transposed of 'jacobian' and 'jacobianInverse', respectively.
      {
        if( (JtAsFieldMatrix - JAsFieldMatrix.transposed()).infinity_norm() != 0 )
        {
          std::cerr << "Error: jacobian and jacobianTransposed are not transposed to each other." << std::endl;
          pass = false;
        }
      }
      {
        if( (JitAsFieldMatrix - JiAsFieldMatrix.transposed()).infinity_norm() != 0 )
        {
          std::cerr << "Error: jacobianInverse and jacobianInverseTransposed are not transposed to each other." << std::endl;
          pass = false;
        }
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
            jtj[ i ][ j ] += JtAsFieldMatrix[ i ][ k ] * JtAsFieldMatrix[ j ][ k ];

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
