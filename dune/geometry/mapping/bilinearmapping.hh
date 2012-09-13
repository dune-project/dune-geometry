// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_BILINEARMAPPING_HH
#define DUNE_GEOMETRY_BILINEARMAPPING_HH

#include <cassert>
#include <limits>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/genericgeometry/geometrytraits.hh>
#include <dune/geometry/genericgeometry/matrixhelper.hh>
#include <dune/geometry/genericgeometry/topologytypes.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class ctype, int dim >
  class ReferenceElement;



  // BilinearMappingTraits
  // ---------------------

  template< class ct >
  struct BilinearMappingTraits
  {
    typedef GenericGeometry::MatrixHelper< GenericGeometry::DuneCoordTraits< ct > > MatrixHelper;

    static ct tolerance () { return 16 * std::numeric_limits< ct >::epsilon(); }

    struct UserData {};
  };



  // BilinearMapping
  // ---------------

  template< class ct, int cdim, class Traits = BilinearMappingTraits< ct > >
  class BilinearMapping
  {
    typedef BilinearMapping< ct, cdim, Traits > This;

  public:
    typedef ct ctype;

    static const int mydimension= 2;
    static const int coorddimension = cdim;

    typedef typename Traits::UserData UserData;

    typedef FieldVector< ctype, mydimension > LocalCoordinate;
    typedef FieldVector< ctype, coorddimension > GlobalCoordinate;

    typedef FieldMatrix< ctype, mydimension, coorddimension > JacobianTransposed;
    typedef FieldMatrix< ctype, coorddimension, mydimension > JacobianInverseTransposed;

    // for compatibility, export the type JacobianInverseTransposed as Jacobian
    typedef JacobianInverseTransposed Jacobian;

    //! type of reference element
    typedef Dune::ReferenceElement< ctype, mydimension > ReferenceElement;

  private:
    typedef typename Traits::MatrixHelper MatrixHelper;

    struct Storage
      : public UserData
    {
      template< class CoordVector >
      Storage ( const CoordVector &coordVector, const UserData &userData )
        : UserData( userData )
      {
        for( int i = 0; i < 4; ++i )
          coefficients[ i ] = coordVector[ i ];
        coefficients[ 1 ] -= coefficients[ 0 ];
        coefficients[ 3 ] -= coefficients[ 2 ];
        coefficients[ 2 ] -= coefficients[ 0 ];
        coefficients[ 3 ] -= coefficients[ 1 ];
      }

      GlobalCoordinate coefficients[ 4 ];
    };

  public:
    template< class CoordVector >
    BilinearMapping ( const CoordVector &coordVector, const UserData &userData = UserData() )
      : storage_( coordVector, userData )
    {}

    template< class CoordVector >
    BilinearMapping ( const ReferenceElement &refElement, const CoordVector &coordVector,
                      const UserData &userData = UserData() )
      : storage_( coordVector, userData )
    {
      assert( refElement.type().isCube() );
    }

    template< class CoordVector >
    BilinearMapping ( Dune::GeometryType gt, const CoordVector &coordVector,
                      const UserData &userData = UserData() )
      : storage_( coordVector, userData )
    {
      assert( (gt.dim() == mydimension) && gt.isCube() );
    }

    /** \brief is this mapping affine? */
    bool affine () const
    {
      return (storage().coefficients[ 3 ].one_norm() < Traits::tolerance());
    }

    /** \brief obtain the name of the reference element */
    Dune::GeometryType type () const
    {
      typedef typename GenericGeometry::CubeTopology< mydimension >::type Topology;
      return Dune::GeometryType( Topology() );
    }

    /** \brief obtain number of corners of the corresponding reference element */
    int numCorners () const { return (1 << mydimension); }

    /** \brief obtain coordinates of the i-th corner */
    GlobalCoordinate corner ( int i ) const
    {
      LocalCoordinate local;
      for( int j = 0; j < mydimension; ++j )
        local[ i ] = ctype( (i >> j) & 1 );
      return global( local );
    }

    /** \brief obtain the centroid of the mapping's image */
    GlobalCoordinate center () const { return global( localBaryCenter() ); }

    /** \brief evaluate the mapping
     *
     *  \param[in]  local  local coordinate to map
     *
     *  \returns corresponding global coordinate
     */
    GlobalCoordinate global ( const LocalCoordinate &local ) const
    {
      GlobalCoordinate global( storage().coefficients[ 0 ] );
      for( int i = 0; i < 2; ++i )
        global.axpy( local[ i ], storage().coefficients[ i+1 ] );
      global.axpy( local[ 0 ]*local[ 1 ], storage().coefficients[ 3 ] );
      return global;
    }

    /** \brief evaluate the inverse mapping
     *
     *  \param[in]  global  global coorindate to map
     *
     *  \return corresponding local coordinate
     *
     *  \note The returned local coordinate y minimizes
     *  \code
     *  (global( x ) - y).two_norm()
     *  \endcode
     */
    LocalCoordinate local ( const GlobalCoordinate &global ) const
    {
      LocalCoordinate local( localBaryCenter() );
      GlobalCoordinate dglobal = (*this).global( local ) - global;

      const ctype tolerance = Traits::tolerance();
      ctype step;
      do {
        step = ctype( 0 );
        for( int j = 0; j < mydimension; ++j )
        {
          GlobalCoordinate direction = derivative( j, local );
          const ctype dd = direction*direction;
          if( dd < tolerance*tolerance )
            continue;

          const ctype alpha = -(dglobal*direction) / dd;
          dglobal.axpy( alpha, direction );
          local[ j ] += alpha;
          step += std::abs( alpha );
        }
      }
      while( step > tolerance );
    }

    /** \brief obtain the integration element
     *
     *  If the Jacobian of the mapping is denoted by $J(x)$, the integration
     *  integration element \f$\mu(x)\f$ is given by
     *  \f[ \mu(x) = \sqrt{|\det (J^T(x) J(x))|}.\f]
     *
     *  \param[in]  local  local coordinate to evaluate the integration element in
     *
     *  \returns the integration element \f$\mu(x)\f$.
     *
     *  \note For affine mappings, it is more efficient to call
     *        jacobianInverseTransposed before integrationElement, if both
     *        are required.
     */
    ctype integrationElement ( const LocalCoordinate &local ) const
    {
      return MatrixHelper::template sqrtDetAAT< mydimension, coorddimension >( jacobianTransosed( local ) );
    }

    /** \brief obtain the volume of the mapping's image
     *
     *  \note The current implementation just returns
     *  \code
     *  integrationElement( baryCenter() ) * ReferenceElement::volume()
     *  \endcode
     *  which is wrong for n-linear surface maps and other nonlinear maps.
     */
    ctype volume () const
    {
      return integrationElement( localBaryCenter() );
    }

    /** \brief obtain the transposed of the Jacobian
     *
     *  \param[in]  local  local coordinate to evaluate Jacobian in
     *
     *  \returns a reference to the transposed of the Jacobian
     *
     *  \note The returned reference is reused on the next call to
     *        JacobianTransposed, destroying the previous value.
     */
    const JacobianTransposed &jacobianTransposed ( const LocalCoordinate &local ) const
    {
      for( int i = 0; i < mydimension; ++i )
        jacobianTransposed_[ i ] = derivative( i, local );
      return jacobianTransposed_;
    }

    /** \brief obtain the transposed of the Jacobian's inverse
     *
     *  The Jacobian's inverse is defined as a pseudo-inverse. If we denote
     *  the Jacobian by \f$J(x)\f$, the following condition holds:
     *  \f[J^{-1}(x) J(x) = I.\f]
     */
    const JacobianInverseTransposed &jacobianInverseTransposed ( const LocalCoordinate &local ) const
    {
      MatrixHelper::template rightInvA< mydimension, coorddimension >( jacobianTransposed(), jacobianInverseTransposed_ );
      return jacobianInverseTransposed_;
    }

    const UserData &userData () const { return storage_; }
    UserData &userData () { return storage_; }

  protected:
    const Storage &storage () const { return storage_; }
    Storage &storage () { return storage_; }

    LocalCoordinate localBaryCenter () const
    {
      LocalCoordinate local;
      for( int i = 0; i < mydimension; ++i )
        local[ i ] = ctype( 1 ) / ctype( 2 );
      return local;
    }

    GlobalCoordinate derivative ( int i, const LocalCoordinate &local ) const
    {
      GlobalCoordinate derivative = storage().coefficients[ i+1 ];
      derivative.axpy( local[ 1-i ], storage().coefficients[ 3 ] );
      return derivative;
    }

  private:
    Storage storage_;
    mutable JacobianTransposed jacobianTransposed_;
    mutable JacobianInverseTransposed jacobianInverseTransposed_;
  };

} // namespace Dune

#endif // #ifndef DUNE_GEOMETRY_BILINEARMAPPING_HH
