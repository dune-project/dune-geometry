// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GEOMETRY_AFFINEGEOMETRY_HH
#define DUNE_GEOMETRY_AFFINEGEOMETRY_HH

/** \file
 *  \brief An implementation of the Geometry interface for affine geometries
 *  \author Martin Nolte
 */

#include <cmath>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  namespace Geo
  {

    template< typename Implementation >
    class ReferenceElement;

    template< class ctype, int dim >
    class ReferenceElementImplementation;

    template< class ctype, int dim >
    struct ReferenceElements;

  }


  namespace Impl
  {

    // FieldMatrixHelper
    // -----------------

    template< class ct >
    struct FieldMatrixHelper
    {
      typedef ct ctype;

      template< int m, int n >
      static void Ax ( const FieldMatrix< ctype, m, n > &A, const FieldVector< ctype, n > &x, FieldVector< ctype, m > &ret )
      {
        for( int i = 0; i < m; ++i )
        {
          ret[ i ] = ctype( 0 );
          for( int j = 0; j < n; ++j )
            ret[ i ] += A[ i ][ j ] * x[ j ];
        }
      }

      template< int m, int n >
      static void ATx ( const FieldMatrix< ctype, m, n > &A, const FieldVector< ctype, m > &x, FieldVector< ctype, n > &ret )
      {
        for( int i = 0; i < n; ++i )
        {
          ret[ i ] = ctype( 0 );
          for( int j = 0; j < m; ++j )
            ret[ i ] += A[ j ][ i ] * x[ j ];
        }
      }

      template< int m, int n, int p >
      static void AB ( const FieldMatrix< ctype, m, n > &A, const FieldMatrix< ctype, n, p > &B, FieldMatrix< ctype, m, p > &ret )
      {
        for( int i = 0; i < m; ++i )
        {
          for( int j = 0; j < p; ++j )
          {
            ret[ i ][ j ] = ctype( 0 );
            for( int k = 0; k < n; ++k )
              ret[ i ][ j ] += A[ i ][ k ] * B[ k ][ j ];
          }
        }
      }

      template< int m, int n, int p >
      static void ATBT ( const FieldMatrix< ctype, m, n > &A, const FieldMatrix< ctype, p, m > &B, FieldMatrix< ctype, n, p > &ret )
      {
        for( int i = 0; i < n; ++i )
        {
          for( int j = 0; j < p; ++j )
          {
            ret[ i ][ j ] = ctype( 0 );
            for( int k = 0; k < m; ++k )
              ret[ i ][ j ] += A[ k ][ i ] * B[ j ][ k ];
          }
        }
      }

      template< int m, int n >
      static void ATA_L ( const FieldMatrix< ctype, m, n > &A, FieldMatrix< ctype, n, n > &ret )
      {
        for( int i = 0; i < n; ++i )
        {
          for( int j = 0; j <= i; ++j )
          {
            ret[ i ][ j ] = ctype( 0 );
            for( int k = 0; k < m; ++k )
              ret[ i ][ j ] += A[ k ][ i ] * A[ k ][ j ];
          }
        }
      }

      template< int m, int n >
      static void ATA ( const FieldMatrix< ctype, m, n > &A, FieldMatrix< ctype, n, n > &ret )
      {
        for( int i = 0; i < n; ++i )
        {
          for( int j = 0; j <= i; ++j )
          {
            ret[ i ][ j ] = ctype( 0 );
            for( int k = 0; k < m; ++k )
              ret[ i ][ j ] += A[ k ][ i ] * A[ k ][ j ];
            ret[ j ][ i ] = ret[ i ][ j ];
          }

          ret[ i ][ i ] = ctype( 0 );
          for( int k = 0; k < m; ++k )
            ret[ i ][ i ] += A[ k ][ i ] * A[ k ][ i ];
        }
      }

      template< int m, int n >
      static void AAT_L ( const FieldMatrix< ctype, m, n > &A, FieldMatrix< ctype, m, m > &ret )
      {
        /*
           if (m==2) {
           ret[0][0] = A[0]*A[0];
           ret[1][1] = A[1]*A[1];
           ret[1][0] = A[0]*A[1];
           }
           else
         */
        for( int i = 0; i < m; ++i )
        {
          for( int j = 0; j <= i; ++j )
          {
            ctype &retij = ret[ i ][ j ];
            retij = A[ i ][ 0 ] * A[ j ][ 0 ];
            for( int k = 1; k < n; ++k )
              retij += A[ i ][ k ] * A[ j ][ k ];
          }
        }
      }

      template< int m, int n >
      static void AAT ( const FieldMatrix< ctype, m, n > &A, FieldMatrix< ctype, m, m > &ret )
      {
        for( int i = 0; i < m; ++i )
        {
          for( int j = 0; j < i; ++j )
          {
            ret[ i ][ j ] = ctype( 0 );
            for( int k = 0; k < n; ++k )
              ret[ i ][ j ] += A[ i ][ k ] * A[ j ][ k ];
            ret[ j ][ i ] = ret[ i ][ j ];
          }
          ret[ i ][ i ] = ctype( 0 );
          for( int k = 0; k < n; ++k )
            ret[ i ][ i ] += A[ i ][ k ] * A[ i ][ k ];
        }
      }

      template< int n >
      static void Lx ( const FieldMatrix< ctype, n, n > &L, const FieldVector< ctype, n > &x, FieldVector< ctype, n > &ret )
      {
        for( int i = 0; i < n; ++i )
        {
          ret[ i ] = ctype( 0 );
          for( int j = 0; j <= i; ++j )
            ret[ i ] += L[ i ][ j ] * x[ j ];
        }
      }

      template< int n >
      static void LTx ( const FieldMatrix< ctype, n, n > &L, const FieldVector< ctype, n > &x, FieldVector< ctype, n > &ret )
      {
        for( int i = 0; i < n; ++i )
        {
          ret[ i ] = ctype( 0 );
          for( int j = i; j < n; ++j )
            ret[ i ] += L[ j ][ i ] * x[ j ];
        }
      }

      template< int n >
      static void LTL ( const FieldMatrix< ctype, n, n > &L, FieldMatrix< ctype, n, n > &ret )
      {
        for( int i = 0; i < n; ++i )
        {
          for( int j = 0; j < i; ++j )
          {
            ret[ i ][ j ] = ctype( 0 );
            for( int k = i; k < n; ++k )
              ret[ i ][ j ] += L[ k ][ i ] * L[ k ][ j ];
            ret[ j ][ i ] = ret[ i ][ j ];
          }
          ret[ i ][ i ] = ctype( 0 );
          for( int k = i; k < n; ++k )
            ret[ i ][ i ] += L[ k ][ i ] * L[ k ][ i ];
        }
      }

      template< int n >
      static void LLT ( const FieldMatrix< ctype, n, n > &L, FieldMatrix< ctype, n, n > &ret )
      {
        for( int i = 0; i < n; ++i )
        {
          for( int j = 0; j < i; ++j )
          {
            ret[ i ][ j ] = ctype( 0 );
            for( int k = 0; k <= j; ++k )
              ret[ i ][ j ] += L[ i ][ k ] * L[ j ][ k ];
            ret[ j ][ i ] = ret[ i ][ j ];
          }
          ret[ i ][ i ] = ctype( 0 );
          for( int k = 0; k <= i; ++k )
            ret[ i ][ i ] += L[ i ][ k ] * L[ i ][ k ];
        }
      }

      template< int n >
      static bool cholesky_L ( const FieldMatrix< ctype, n, n > &A, FieldMatrix< ctype, n, n > &ret, const bool checkSingular = false )
      {
        using std::sqrt;
        for( int i = 0; i < n; ++i )
        {
          ctype &rii = ret[ i ][ i ];

          ctype xDiag = A[ i ][ i ];
          for( int j = 0; j < i; ++j )
            xDiag -= ret[ i ][ j ] * ret[ i ][ j ];

          // in some cases A can be singular, e.g. when checking local for
          // outside points during checkInside
          if( checkSingular && ! ( xDiag > ctype( 0 )) )
            return false ;

          // otherwise this should be true always
          assert( xDiag > ctype( 0 ) );
          rii = sqrt( xDiag );

          ctype invrii = ctype( 1 ) / rii;
          for( int k = i+1; k < n; ++k )
          {
            ctype x = A[ k ][ i ];
            for( int j = 0; j < i; ++j )
              x -= ret[ i ][ j ] * ret[ k ][ j ];
            ret[ k ][ i ] = invrii * x;
          }
        }

        // return true for meaning A is non-singular
        return true;
      }

      template< int n >
      static ctype detL ( const FieldMatrix< ctype, n, n > &L )
      {
        ctype det( 1 );
        for( int i = 0; i < n; ++i )
          det *= L[ i ][ i ];
        return det;
      }

      template< int n >
      static ctype invL ( FieldMatrix< ctype, n, n > &L )
      {
        ctype det( 1 );
        for( int i = 0; i < n; ++i )
        {
          ctype &lii = L[ i ][ i ];
          det *= lii;
          lii = ctype( 1 ) / lii;
          for( int j = 0; j < i; ++j )
          {
            ctype &lij = L[ i ][ j ];
            ctype x = lij * L[ j ][ j ];
            for( int k = j+1; k < i; ++k )
              x += L[ i ][ k ] * L[ k ][ j ];
            lij = (-lii) * x;
          }
        }
        return det;
      }

      // calculates x := L^{-1} x
      template< int n >
      static void invLx ( FieldMatrix< ctype, n, n > &L, FieldVector< ctype, n > &x )
      {
        for( int i = 0; i < n; ++i )
        {
          for( int j = 0; j < i; ++j )
            x[ i ] -= L[ i ][ j ] * x[ j ];
          x[ i ] /= L[ i ][ i ];
        }
      }

      // calculates x := L^{-T} x
      template< int n >
      static void invLTx ( FieldMatrix< ctype, n, n > &L, FieldVector< ctype, n > &x )
      {
        for( int i = n; i > 0; --i )
        {
          for( int j = i; j < n; ++j )
            x[ i-1 ] -= L[ j ][ i-1 ] * x[ j ];
          x[ i-1 ] /= L[ i-1 ][ i-1 ];
        }
      }

      template< int n >
      static ctype spdDetA ( const FieldMatrix< ctype, n, n > &A )
      {
        // return A[0][0]*A[1][1]-A[1][0]*A[1][0];
        FieldMatrix< ctype, n, n > L;
        cholesky_L( A, L );
        return detL( L );
      }

      template< int n >
      static ctype spdInvA ( FieldMatrix< ctype, n, n > &A )
      {
        FieldMatrix< ctype, n, n > L;
        cholesky_L( A, L );
        const ctype det = invL( L );
        LTL( L, A );
        return det;
      }

      // calculate x := A^{-1} x
      template< int n >
      static bool spdInvAx ( FieldMatrix< ctype, n, n > &A, FieldVector< ctype, n > &x, const bool checkSingular = false )
      {
        FieldMatrix< ctype, n, n > L;
        const bool invertible = cholesky_L( A, L, checkSingular );
        if( ! invertible ) return invertible ;
        invLx( L, x );
        invLTx( L, x );
        return invertible;
      }

      template< int m, int n >
      static ctype detATA ( const FieldMatrix< ctype, m, n > &A )
      {
        if( m >= n )
        {
          FieldMatrix< ctype, n, n > ata;
          ATA_L( A, ata );
          return spdDetA( ata );
        }
        else
          return ctype( 0 );
      }

      /** \brief Compute the square root of the determinant of A times A transposed
       *
       *  This is the volume element for an embedded submanifold and needed to
       *  implement the method integrationElement().
       */
      template< int m, int n >
      static ctype sqrtDetAAT ( const FieldMatrix< ctype, m, n > &A )
      {
        using std::abs;
        using std::sqrt;
        // These special cases are here not only for speed reasons:
        // The general implementation aborts if the matrix is almost singular,
        // and the special implementation provide a stable way to handle that case.
        if( (n == 2) && (m == 2) )
        {
          // Special implementation for 2x2 matrices: faster and more stable
          return abs( A[ 0 ][ 0 ]*A[ 1 ][ 1 ] - A[ 1 ][ 0 ]*A[ 0 ][ 1 ] );
        }
        else if( (n == 3) && (m == 3) )
        {
          // Special implementation for 3x3 matrices
          const ctype v0 = A[ 0 ][ 1 ] * A[ 1 ][ 2 ] - A[ 1 ][ 1 ] * A[ 0 ][ 2 ];
          const ctype v1 = A[ 0 ][ 2 ] * A[ 1 ][ 0 ] - A[ 1 ][ 2 ] * A[ 0 ][ 0 ];
          const ctype v2 = A[ 0 ][ 0 ] * A[ 1 ][ 1 ] - A[ 1 ][ 0 ] * A[ 0 ][ 1 ];
          return abs( v0 * A[ 2 ][ 0 ] + v1 * A[ 2 ][ 1 ] + v2 * A[ 2 ][ 2 ] );
        }
        else if ( (n == 3) && (m == 2) )
        {
          // Special implementation for 2x3 matrices
          const ctype v0 = A[ 0 ][ 0 ] * A[ 1 ][ 1 ] - A[ 0 ][ 1 ] * A[ 1 ][ 0 ];
          const ctype v1 = A[ 0 ][ 0 ] * A[ 1 ][ 2 ] - A[ 1 ][ 0 ] * A[ 0 ][ 2 ];
          const ctype v2 = A[ 0 ][ 1 ] * A[ 1 ][ 2 ] - A[ 0 ][ 2 ] * A[ 1 ][ 1 ];
          return sqrt( v0*v0 + v1*v1 + v2*v2);
        }
        else if( n >= m )
        {
          // General case
          FieldMatrix< ctype, m, m > aat;
          AAT_L( A, aat );
          return spdDetA( aat );
        }
        else
          return ctype( 0 );
      }

      // A^{-1}_L = (A^T A)^{-1} A^T
      // => A^{-1}_L A = I
      template< int m, int n >
      static ctype leftInvA ( const FieldMatrix< ctype, m, n > &A, FieldMatrix< ctype, n, m > &ret )
      {
        static_assert((m >= n), "Matrix has no left inverse.");
        FieldMatrix< ctype, n, n > ata;
        ATA_L( A, ata );
        const ctype det = spdInvA( ata );
        ATBT( ata, A, ret );
        return det;
      }

      template< int m, int n >
      static void leftInvAx ( const FieldMatrix< ctype, m, n > &A, const FieldVector< ctype, m > &x, FieldVector< ctype, n > &y )
      {
        static_assert((m >= n), "Matrix has no left inverse.");
        FieldMatrix< ctype, n, n > ata;
        ATx( A, x, y );
        ATA_L( A, ata );
        spdInvAx( ata, y );
      }

      /** \brief Compute right pseudo-inverse of matrix A */
      template< int m, int n >
      static ctype rightInvA ( const FieldMatrix< ctype, m, n > &A, FieldMatrix< ctype, n, m > &ret )
      {
        static_assert((n >= m), "Matrix has no right inverse.");
        using std::abs;
        if( (n == 2) && (m == 2) )
        {
          const ctype det = (A[ 0 ][ 0 ]*A[ 1 ][ 1 ] - A[ 1 ][ 0 ]*A[ 0 ][ 1 ]);
          const ctype detInv = ctype( 1 ) / det;
          ret[ 0 ][ 0 ] = A[ 1 ][ 1 ] * detInv;
          ret[ 1 ][ 1 ] = A[ 0 ][ 0 ] * detInv;
          ret[ 1 ][ 0 ] = -A[ 1 ][ 0 ] * detInv;
          ret[ 0 ][ 1 ] = -A[ 0 ][ 1 ] * detInv;
          return abs( det );
        }
        else
        {
          FieldMatrix< ctype, m , m > aat;
          AAT_L( A, aat );
          const ctype det = spdInvA( aat );
          ATBT( A , aat , ret );
          return det;
        }
      }

      template< int m, int n >
      static bool xTRightInvA ( const FieldMatrix< ctype, m, n > &A, const FieldVector< ctype, n > &x, FieldVector< ctype, m > &y )
      {
        static_assert((n >= m), "Matrix has no right inverse.");
        FieldMatrix< ctype, m, m > aat;
        Ax( A, x, y );
        AAT_L( A, aat );
        // check whether aat is singular and return true if non-singular
        return spdInvAx( aat, y, true );
      }
    };

  } // namespace Impl



  /** \brief Implementation of the Geometry interface for affine geometries
   * \tparam ct Type used for coordinates
   * \tparam mydim Dimension of the geometry
   * \tparam cdim Dimension of the world space
   */
  template< class ct, int mydim, int cdim>
  class AffineGeometry
  {
  public:

    /** \brief Type used for coordinates */
    typedef ct ctype;

    /** \brief Dimension of the geometry */
    static const int mydimension= mydim;

    /** \brief Dimension of the world space */
    static const int coorddimension = cdim;

    /** \brief Type for local coordinate vector */
    typedef FieldVector< ctype, mydimension > LocalCoordinate;

    /** \brief Type for coordinate vector in world space */
    typedef FieldVector< ctype, coorddimension > GlobalCoordinate;

    /** \brief Type used for volume */
    typedef ctype Volume;

    /** \brief Type for the transposed Jacobian matrix */
    typedef FieldMatrix< ctype, mydimension, coorddimension > JacobianTransposed;

    /** \brief Type for the transposed inverse Jacobian matrix */
    typedef FieldMatrix< ctype, coorddimension, mydimension > JacobianInverseTransposed;

    /** \brief Type for the Jacobian matrix */
    typedef FieldMatrix< ctype, coorddimension, mydimension > Jacobian;

    /** \brief Type for the inverse Jacobian matrix */
    typedef FieldMatrix< ctype, mydimension, coorddimension > JacobianInverse;

  private:
    //! type of reference element
    typedef Geo::ReferenceElement< Geo::ReferenceElementImplementation< ctype, mydimension > > ReferenceElement;

    typedef Geo::ReferenceElements< ctype, mydimension > ReferenceElements;

    // Helper class to compute a matrix pseudo inverse
    typedef Impl::FieldMatrixHelper< ct > MatrixHelper;

  public:
    /** \brief Create affine geometry from reference element, one vertex, and the Jacobian matrix */
    AffineGeometry ( const ReferenceElement &refElement, const GlobalCoordinate &origin,
                     const JacobianTransposed &jt )
      : refElement_(refElement), origin_(origin), jacobianTransposed_(jt)
    {
      integrationElement_ = MatrixHelper::template rightInvA< mydimension, coorddimension >( jacobianTransposed_, jacobianInverseTransposed_ );
    }

    /** \brief Create affine geometry from GeometryType, one vertex, and the Jacobian matrix */
    AffineGeometry ( Dune::GeometryType gt, const GlobalCoordinate &origin,
                     const JacobianTransposed &jt )
      : AffineGeometry(ReferenceElements::general( gt ), origin, jt)
    { }

    /** \brief Create affine geometry from reference element and a vector of vertex coordinates */
    template< class CoordVector >
    AffineGeometry ( const ReferenceElement &refElement, const CoordVector &coordVector )
      : refElement_(refElement), origin_(coordVector[0])
    {
      for( int i = 0; i < mydimension; ++i )
        jacobianTransposed_[ i ] = coordVector[ i+1 ] - origin_;
      integrationElement_ = MatrixHelper::template rightInvA< mydimension, coorddimension >( jacobianTransposed_, jacobianInverseTransposed_ );
    }

    /** \brief Create affine geometry from GeometryType and a vector of vertex coordinates */
    template< class CoordVector >
    AffineGeometry ( Dune::GeometryType gt, const CoordVector &coordVector )
      : AffineGeometry(ReferenceElements::general( gt ), coordVector)
    { }

    /** \brief Always true: this is an affine geometry */
    bool affine () const { return true; }

    /** \brief Obtain the type of the reference element */
    Dune::GeometryType type () const { return refElement_.type(); }

    /** \brief Obtain number of corners of the corresponding reference element */
    int corners () const { return refElement_.size( mydimension ); }

    /** \brief Obtain coordinates of the i-th corner */
    GlobalCoordinate corner ( int i ) const
    {
      return global( refElement_.position( i, mydimension ) );
    }

    /** \brief Obtain the centroid of the mapping's image */
    GlobalCoordinate center () const { return global( refElement_.position( 0, 0 ) ); }

    /** \brief Evaluate the mapping
     *
     *  \param[in]  local  local coordinate to map
     *
     *  \returns corresponding global coordinate
     */
    GlobalCoordinate global ( const LocalCoordinate &local ) const
    {
      GlobalCoordinate global( origin_ );
      jacobianTransposed_.umtv( local, global );
      return global;
    }

    /** \brief Evaluate the inverse mapping
     *
     *  \param[in]  global  global coordinate to map
     *
     *  \return corresponding local coordinate
     *
     *  The returned local coordinate y minimizes
     *  \code
     *  (global( y ) - x).two_norm()
     *  \endcode
     *  on the entire affine hull of the reference element.  This degenerates
     *  to the inverse map if the argument y is in the range of the map.
     */
    LocalCoordinate local ( const GlobalCoordinate &global ) const
    {
      LocalCoordinate local;
      jacobianInverseTransposed_.mtv( global - origin_, local );
      return local;
    }

    /** \brief Obtain the integration element
     *
     *  If the Jacobian of the mapping is denoted by $J(x)$, the integration
     *  integration element \f$\mu(x)\f$ is given by
     *  \f[ \mu(x) = \sqrt{|\det (J^T(x) J(x))|}.\f]
     *
     *  \param[in]  local  local coordinate to evaluate the integration element in
     *
     *  \returns the integration element \f$\mu(x)\f$.
     */
    ctype integrationElement ([[maybe_unused]] const LocalCoordinate &local) const
    {
      return integrationElement_;
    }

    /** \brief Obtain the volume of the element */
    Volume volume () const
    {
      return integrationElement_ * refElement_.volume();
    }

    /** \brief Obtain the transposed of the Jacobian
     *
     *  \param[in]  local  local coordinate to evaluate Jacobian in
     *
     *  \returns a reference to the transposed of the Jacobian
     */
    const JacobianTransposed &jacobianTransposed ([[maybe_unused]] const LocalCoordinate &local) const
    {
      return jacobianTransposed_;
    }

    /** \brief Obtain the transposed of the Jacobian's inverse
     *
     *  The Jacobian's inverse is defined as a pseudo-inverse. If we denote
     *  the Jacobian by \f$J(x)\f$, the following condition holds:
     *  \f[J^{-1}(x) J(x) = I.\f]
     */
    const JacobianInverseTransposed &jacobianInverseTransposed ([[maybe_unused]] const LocalCoordinate &local) const
    {
      return jacobianInverseTransposed_;
    }

    /** \brief Obtain the Jacobian
     *
     *  \param[in]  local  local coordinate to evaluate Jacobian in
     *
     *  \returns a copy of the transposed of the Jacobian
     */
    Jacobian jacobian ([[maybe_unused]] const LocalCoordinate &local) const
    {
      return jacobianTransposed_.transposed();
    }

    /** \brief Obtain the Jacobian's inverse
     *
     *  The Jacobian's inverse is defined as a pseudo-inverse. If we denote
     *  the Jacobian by \f$J(x)\f$, the following condition holds:
     *  \f[J^{-1}(x) J(x) = I.\f]
     */
    JacobianInverse jacobianInverse ([[maybe_unused]] const LocalCoordinate &local) const
    {
      return jacobianInverseTransposed_.transposed();
    }

    friend ReferenceElement referenceElement ( const AffineGeometry &geometry )
    {
      return geometry.refElement_;
    }

  private:
    ReferenceElement refElement_;
    GlobalCoordinate origin_;
    JacobianTransposed jacobianTransposed_;
    JacobianInverseTransposed jacobianInverseTransposed_;
    ctype integrationElement_;
  };

} // namespace Dune

#endif // #ifndef DUNE_GEOMETRY_AFFINEGEOMETRY_HH
