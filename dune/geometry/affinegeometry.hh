// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_AFFINEGEOMETRY_HH
#define DUNE_GEOMETRY_AFFINEGEOMETRY_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/genericgeometry/geometrytraits.hh>
#include <dune/geometry/genericgeometry/matrixhelper.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class ctype, int dim >
  class ReferenceElement;

  template< class ctype, int dim >
  class ReferenceElements;



  // AffineGeometryTraits
  // --------------------

  template< class ct >
  struct AffineGeometryTraits
  {
    typedef GenericGeometry::MatrixHelper< GenericGeometry::DuneCoordTraits< ct > > MatrixHelper;
  };



  // AffineGeometry
  // --------------

  template< class ct, int mydim, int cdim, class Traits = AffineGeometryTraits< ct > >
  class AffineGeometry
  {
    typedef AffineGeometry< ct, mydim, cdim, Traits > This;

  public:
    typedef ct ctype;

    static const int mydimension= mydim;
    static const int coorddimension = cdim;

    typedef FieldVector< ctype, mydimension > LocalCoordinate;
    typedef FieldVector< ctype, coorddimension > GlobalCoordinate;

    typedef FieldMatrix< ctype, mydimension, coorddimension > JacobianTransposed;
    typedef FieldMatrix< ctype, coorddimension, mydimension > JacobianInverseTransposed;

    // for compatibility, export the type JacobianInverseTransposed as Jacobian
    typedef JacobianInverseTransposed Jacobian;

    //! type of reference element
    typedef Dune::ReferenceElement< ctype, mydimension > ReferenceElement;

  protected:
    typedef Dune::ReferenceElements< ctype, mydimension > ReferenceElements;

  private:
    typedef typename Traits::MatrixHelper MatrixHelper;

    struct Storage
    {
      Storage ( const ReferenceElement &refEl, const GlobalCoordinate &org,
                const JacobianTransposed &jt )
        : refElement( &refEl ),
          origin( org ),
          jacobianTransposed( jt )
      {
        integrationElement = MatrixHelper::template rightInvA< mydimension, coorddimension >( jacobianTransposed, jacobianInverseTransposed );
      }

      template< class CoordVector >
      Storage ( const ReferenceElement &refEl, const CoordVector &coordVector)
        : refElement( &refEl ),
          origin( coordVector[ 0 ] )
      {
        for( int i = 0; i < mydimension; ++i )
          jacobianTransposed[ i ] = coordVector[ i+1 ] - origin;
        integrationElement = MatrixHelper::template rightInvA< mydimension, coorddimension >( jacobianTransposed, jacobianInverseTransposed );
      }

      const ReferenceElement *refElement;
      GlobalCoordinate origin;
      JacobianTransposed jacobianTransposed;
      JacobianInverseTransposed jacobianInverseTransposed;
      ctype integrationElement;
    };

  public:
    AffineGeometry ( const ReferenceElement &refElement, const GlobalCoordinate &origin,
                     const JacobianTransposed &jt )
      : storage_( refElement, origin, jt )
    {}

    AffineGeometry ( Dune::GeometryType gt, const GlobalCoordinate &origin,
                     const JacobianTransposed &jt )
      : storage_( ReferenceElements::general( gt ), origin, jt )
    {}

    template< class CoordVector >
    AffineGeometry ( const ReferenceElement &refElement, const CoordVector &coordVector )
      : storage_( refElement, coordVector )
    {}

    template< class CoordVector >
    AffineGeometry ( Dune::GeometryType gt, const CoordVector &coordVector )
      : storage_( ReferenceElements::general( gt ), coordVector )
    {}

    /** \brief is this mapping affine? */
    bool affine () const { return true; }

    /** \brief obtain the name of the reference element */
    Dune::GeometryType type () const { return refElement().type(); }

    /** \brief obtain number of corners of the corresponding reference element */
    int corners () const { return refElement().size( mydimension ); }

    /** \brief obtain coordinates of the i-th corner */
    GlobalCoordinate corner ( int i ) const
    {
      return global( refElement().position( i, mydimension ) );
    }

    /** \brief obtain the centroid of the mapping's image */
    GlobalCoordinate center () const { return global( refElement().position( 0, 0 ) ); }

    /** \brief evaluate the mapping
     *
     *  \param[in]  local  local coordinate to map
     *
     *  \returns corresponding global coordinate
     */
    GlobalCoordinate global ( const LocalCoordinate &local ) const
    {
      GlobalCoordinate global( storage().origin );
      storage().jacobianTransposed.umtv( local, global );
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
      LocalCoordinate local;
      storage().jacobianInverseTransposed.mtv( global - storage().origin, local );
      return local;
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
     */
    ctype integrationElement ( const LocalCoordinate &local ) const
    {
      return storage().integrationElement;
    }

    /** \brief obtain the volume of the mapping's image */
    ctype volume () const
    {
      return storage().integrationElement * refElement().volume();
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
      return storage().jacobianTransposed;
    }

    /** \brief obtain the transposed of the Jacobian's inverse
     *
     *  The Jacobian's inverse is defined as a pseudo-inverse. If we denote
     *  the Jacobian by \f$J(x)\f$, the following condition holds:
     *  \f[J^{-1}(x) J(x) = I.\f]
     */
    const JacobianInverseTransposed &jacobianInverseTransposed ( const LocalCoordinate &local ) const
    {
      return storage().jacobianInverseTransposed;
    }

  protected:
    const Storage &storage () const { return storage_; }
    Storage &storage () { return storage_; }

    const ReferenceElement &refElement () const { return *storage().refElement; }

  private:
    Storage storage_;
  };

} // namespace Dune

#endif // #ifndef DUNE_GEOMETRY_AFFINEGEOMETRY_HH
